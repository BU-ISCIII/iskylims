#!/usr/bin/bash

ISKYLIMS_VERSION="3.1.0"

usage() {
cat << EOF
This script installs and upgrades the iskylims app.

Usage : $0 [--demo_data] [--git_revision] [--compose_file] [--install_conf] [--action] [--script] [--script_before] [--script_after] [--engine] [--test]
    Optional input data:
    --demo_data         | Provide already downloaded demo data from Zenodo
    --git_revision      | Specify the Git revision to install (default: main, or 'current' to use copied local sources)
    --compose_file      | Compose file to use (overrides default)
    --install_conf      | Settings file consumed during container image build (mandatory for production)
    --install_conf_map  | Service-specific settings file: service,path (can be repeated)
    --action            | install (default), upgrade, or fix-permissions
    --script            | Run a Django migration script after migrations (can be repeated)
    --script_before     | Run a Django migration script before migrations (can be repeated)
    --script_after      | Run a Django migration script after migrations (can be repeated)
    --skip_demo_data    | Skip downloading/copying demo data to samba container
    --skip_test_data    | Skip loading test fixtures (test/test_data.json)
    --engine            | Container engine to use: docker (default) or podman
    --test              | Use development/test compose file and sample data

Examples:
    Deploy production container pointing to an external DB/Samba:
    bash $0 --install_conf conf/my_prod_settings_iskylims.txt

    Deploy production with service-specific settings mapping:
    bash $0 --install_conf_map app,conf/docker_production_settings.txt

    Upgrade an existing production deployment using the same database:
    bash $0 --install_conf conf/my_prod_settings_iskylims.txt --action upgrade

    Repair production bind mount and volume permissions without rebuilding or bootstrapping:
    bash $0 --install_conf conf/my_prod_settings_iskylims.txt --action fix-permissions

    Install demo container system with local services
    bash $0 --test

    Install test stack from current local committed sources without checking out a branch in-container
    bash $0 --test --git_revision current

    Provide already downloaded data from Zenodo (compressed) for test environment
    bash $0 --demo_data /path/to/iskylims_demo_data.tar.gz

EOF
}

# translate long options to short
reset=true

for arg in "$@"
do
    if [ -n "$reset" ]; then
      unset reset
      set --      # this resets the "$@" array so we can rebuild it
    fi
    case "$arg" in
        # OPTIONAL
        --demo_data)         set -- "$@" -d ;;
        --git_revision)      set -- "$@" -g ;;
        --compose_file)      set -- "$@" -c ;;
        --install_conf)      set -- "$@" -s ;;
        --install_conf_map)  set -- "$@" -j ;;
        --action)            set -- "$@" -a ;;
        --script)            set -- "$@" -m ;;
        --script_before)     set -- "$@" -b ;;
        --script_after)      set -- "$@" -f ;;
        --skip_demo_data)    set -- "$@" -n ;;
        --skip_test_data)    set -- "$@" -t ;;
        --test)              set -- "$@" -p ;;
        --engine)            set -- "$@" -e ;;

        # ADDITIONAL
        --help)              set -- "$@" -h ;;
        --version)           set -- "$@" -v ;;
        # PASSING VALUE IN PARAMETER
        *)                   set -- "$@" "$arg" ;;
    esac
done

# SETTING DEFAULT VALUES
demo_data=false
git_revision="main"
compose_file=""
install_conf=""
install_conf_container=""
install_conf_map_entries=()
skip_demo_data=""
skip_test_data=""
mode="production"
action="install"
run_script=false
run_script_before=false
migration_script=()
migration_script_before=()
engine="docker"

ENGINE_CMD=()
COMPOSE_CMD=()

set_engine() {
    if [ "$engine" = "docker" ]; then
        if ! command -v docker >/dev/null 2>&1; then
            echo "docker not found. Install docker or use --engine podman."
            exit 1
        fi
        ENGINE_CMD=("docker")
        COMPOSE_CMD=("docker" "compose")
    else
        if ! command -v podman >/dev/null 2>&1; then
            echo "podman not found. Install podman or use --engine docker."
            exit 1
        fi
        ENGINE_CMD=("podman")
        if command -v podman-compose >/dev/null 2>&1; then
            COMPOSE_CMD=("podman-compose")
        elif podman compose version >/dev/null 2>&1; then
            COMPOSE_CMD=("podman" "compose")
        else
            echo "podman compose not available. Install podman-compose or use --engine docker."
            exit 1
        fi
    fi
}

engine_exec() {
    "${ENGINE_CMD[@]}" "$@"
}

compose_exec() {
    "${COMPOSE_CMD[@]}" "$@"
}

copy_with_podman_fallback() {
    local src="$1"
    local dst="$2"

    if cp "$src" "$dst" 2>/dev/null; then
        return 0
    fi

    if [ "$engine" = "podman" ]; then
        if podman unshare cp "$src" "$dst"; then
            return 0
        fi
    fi

    echo "Failed to copy '$src' to '$dst'" >&2
    return 1
}

chmod_with_podman_fallback() {
    local mode="$1"
    shift

    if chmod "$mode" "$@" 2>/dev/null; then
        return 0
    fi

    if [ "$engine" = "podman" ]; then
        if podman unshare chmod "$mode" "$@"; then
            return 0
        fi
    fi

    echo "Failed to chmod $mode: $*" >&2
    return 1
}

chown_with_podman_fallback() {
    local owner="$1"
    shift

    if chown -R "$owner" "$@" 2>/dev/null; then
        return 0
    fi

    if [ "$engine" = "podman" ]; then
        if podman unshare chown -R "$owner" "$@"; then
            return 0
        fi
    fi

    echo "Failed to chown $owner: $*" >&2
    return 1
}

normalize_apache_server_name() {
    local value="$1"

    value="${value#http://}"
    value="${value#https://}"
    value="${value%%/*}"
    value="${value%%:*}"

    if [ -z "$value" ] || [ "$value" = "*" ]; then
        value="localhost"
    fi

    echo "$value"
}

generate_django_secret_key() {
    if command -v python3 >/dev/null 2>&1; then
        python3 -c "import secrets; print(''.join(secrets.choice('abcdefghijklmnopqrstuvwxyz0123456789!@#$%^&*(-_=+)') for _ in range(50)))"
    else
        LC_ALL=C tr -dc 'A-Za-z0-9!@#$%^&*(-_=+)' < /dev/urandom | head -c 50
        printf "\n"
    fi
}

sed_replacement_escape() {
    printf '%s' "$1" | sed -e 's/[\\&|]/\\&/g'
}

render_django_settings_file() {
    local settings_path="$1"
    local secret_line=""
    local tmp_file=""
    local db_user db_pass db_name db_host db_port
    local email_host email_port email_user email_pass email_tls
    local local_server_ip dns_url

    if [ -f "$settings_path" ]; then
        secret_line="$(grep -E "^SECRET_KEY[[:space:]]*=" "$settings_path" | tail -n 1)"
    fi
    if [ -z "$secret_line" ] || [[ "$secret_line" =~ SECRET_KEY[[:space:]]*=[[:space:]]*SECRET ]]; then
        secret_line="SECRET_KEY = '$(generate_django_secret_key)'"
    fi

    db_user="$(read_install_conf_value DB_USER "$host_install_conf_path")"
    db_pass="$(read_install_conf_value DB_PASS "$host_install_conf_path")"
    db_name="$(read_install_conf_value DB_NAME "$host_install_conf_path")"
    db_host="$(read_install_conf_value DB_SERVER_IP "$host_install_conf_path")"
    db_port="$(read_install_conf_value DB_PORT "$host_install_conf_path")"
    email_host="$(read_install_conf_value EMAIL_HOST_SERVER "$host_install_conf_path")"
    email_port="$(read_install_conf_value EMAIL_PORT "$host_install_conf_path")"
    email_user="$(read_install_conf_value EMAIL_HOST_USER "$host_install_conf_path")"
    email_pass="$(read_install_conf_value EMAIL_HOST_PASSWORD "$host_install_conf_path")"
    email_tls="$(read_install_conf_value EMAIL_USE_TLS "$host_install_conf_path")"
    local_server_ip="$(read_install_conf_value LOCAL_SERVER_IP "$host_install_conf_path")"
    dns_url="$(read_install_conf_value DNS_URL "$host_install_conf_path")"

    tmp_file="$(mktemp)"
    cp "$repo_root/conf/template_settings.txt" "$tmp_file"
    sed -i \
        -e "s|^SECRET_KEY.*|$(sed_replacement_escape "$secret_line")|" \
        -e "s|djangouser|$(sed_replacement_escape "$db_user")|g" \
        -e "s|djangopass|$(sed_replacement_escape "$db_pass")|g" \
        -e "s|djangohost|$(sed_replacement_escape "$db_host")|g" \
        -e "s|djangoport|$(sed_replacement_escape "$db_port")|g" \
        -e "s|djangodbname|$(sed_replacement_escape "$db_name")|g" \
        -e "s|emailhostserver|$(sed_replacement_escape "$email_host")|g" \
        -e "s|emailport|$(sed_replacement_escape "$email_port")|g" \
        -e "s|emailhostuser|$(sed_replacement_escape "$email_user")|g" \
        -e "s|emailhostpassword|$(sed_replacement_escape "$email_pass")|g" \
        -e "s|emailhosttls|$(sed_replacement_escape "$email_tls")|g" \
        -e "s|localserverip|$(sed_replacement_escape "$local_server_ip")|g" \
        -e "s|localhost|$(sed_replacement_escape "$dns_url")|g" \
        "$tmp_file"

    if copy_with_podman_fallback "$tmp_file" "$settings_path"; then
        if ! chmod_with_podman_fallback 0664 "$settings_path"; then
            rm -f "$tmp_file"
            return 1
        fi
        rm -f "$tmp_file"
        return 0
    fi

    rm -f "$tmp_file"
    return 1
}

normalize_settings_bind_path() {
    local value="$1"

    if [ -z "$value" ]; then
        echo "$install_path/iskylims/settings.py"
        return 0
    fi

    if [ -d "$value" ] || [[ "$value" = */ ]] || [[ "$value" != *.py ]]; then
        echo "${value%/}/settings.py"
        return 0
    fi

    echo "$value"
}

render_apache_config() {
    local src="$1"
    local dst="$2"
    local tmp_file=""

    tmp_file="$(mktemp)"
    sed \
        -e "s|__ISKYLIMS_SERVER_NAME__|$apache_server_name|g" \
        -e "s|__ISKYLIMS_LOG_NAME__|$apache_log_name|g" \
        -e "s|__INSTALL_PATH__|$install_path|g" \
        -e "s|__APP_PORT__|$app_port|g" \
        -e "s|__ISKYLIMS_FORWARDED_PROTO__|$(sed_replacement_escape "$apache_forwarded_proto")|g" \
        -e "s|__ISKYLIMS_FORWARDED_PORT__|$(sed_replacement_escape "$apache_forwarded_port")|g" \
        -e "s|__GUNICORN_TIMEOUT__|$gunicorn_timeout|g" \
        -e "s|__SERVER_STATUS_SERVER_NAME__|$(sed_replacement_escape "$apache_status_server_name")|g" \
        -e "s|__SERVER_STATUS_ALIASES__|$(sed_replacement_escape "$apache_status_aliases")|g" \
        -e "s|__SERVER_STATUS_ALLOW_FROM__|$(sed_replacement_escape "$apache_status_allow_from")|g" \
        "$src" > "$tmp_file"

    if copy_with_podman_fallback "$tmp_file" "$dst"; then
        if ! chmod_with_podman_fallback 0664 "$dst"; then
            rm -f "$tmp_file"
            return 1
        fi
        rm -f "$tmp_file"
        return 0
    fi

    rm -f "$tmp_file"
    return 1
}

prepare_django_settings_bind_mount() {
    local settings_path="$1"

    if [ "$mode" != "production" ]; then
        return 0
    fi

    if [ -d "$settings_path" ]; then
        echo "DJANGO_SETTINGS_PATH must resolve to a file path, but '$settings_path' is a directory." >&2
        echo "Use a full path like '$settings_path/settings.py' or remove the directory and rerun." >&2
        return 1
    fi

    mkdir -p "$(dirname "$settings_path")"
    if [ ! -f "$settings_path" ] || grep -Eq "SECRET_KEY[[:space:]]*=[[:space:]]*SECRET|emailhosttls|djangouser|djangopass|djangohost|djangodbname" "$settings_path"; then
        render_django_settings_file "$settings_path"
    fi
    chmod_with_podman_fallback 0664 "$settings_path"
}

# PARSE VARIABLE ARGUMENTS WITH getopts
options=":d:g:c:s:j:a:m:b:f:e:vhntp"
while getopts $options opt; do
    case $opt in
        d)
            demo_data=$OPTARG
            ;;
        g)
            git_revision=$OPTARG
            ;;
        c)
            compose_file=$OPTARG
            ;;
        s)
            install_conf=$OPTARG
            ;;
        j)
            install_conf_map_entries+=("$OPTARG")
            ;;
        a)
            action=$OPTARG
            if [[ "$action" != "install" && "$action" != "upgrade" && "$action" != "fix-permissions" ]]; then
                echo "Invalid action '$action'. Use install, upgrade, or fix-permissions."
                exit 1
            fi
            ;;
        m)
            run_script=true
            migration_script+=("$OPTARG")
            ;;
        b)
            run_script_before=true
            migration_script_before+=("$OPTARG")
            ;;
        e)
            engine=$OPTARG
            if [[ "$engine" != "docker" && "$engine" != "podman" ]]; then
                echo "Invalid engine '$engine'. Use docker or podman."
                exit 1
            fi
            ;;
        f)
            run_script=true
            migration_script+=("$OPTARG")
            ;;
        n)
            skip_demo_data=true
            ;;
        t)
            skip_test_data=true
            ;;
        p)
            mode="test"
            ;;
        h)
            usage
            exit 1
            ;;
        v)
            echo $ISKYLIMS_VERSION
            exit 1
            ;;
        \?)
            echo "Invalid Option: -$OPTARG" 1>&2
            usage
            exit 1
            ;;
        : )
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
        * )
            echo "Unimplemented option: -$OPTARG" >&2;
            exit 1
            ;;
    esac
done
shift $((OPTIND-1))

if [ "$mode" = "test" ]; then
    if [ -z "$compose_file" ]; then
        compose_file="docker-compose.test.yml"
    fi
else
    if [ -z "$compose_file" ]; then
        compose_file="docker-compose.prod.yml"
    fi
fi

app_service="${APP_SERVICE:-app}"
selected_install_conf="$install_conf"
for map_entry in "${install_conf_map_entries[@]}"; do
    svc_name="${map_entry%%,*}"
    conf_name="${map_entry#*,}"
    if [ -z "$svc_name" ] || [ -z "$conf_name" ] || [ "$svc_name" = "$map_entry" ]; then
        echo "Invalid --install_conf_map value '$map_entry'. Expected format: service,path"
        exit 1
    fi
    if [ "$svc_name" != "app" ]; then
        echo "Unknown service '$svc_name' in --install_conf_map. Valid service: app"
        exit 1
    fi
    selected_install_conf="$conf_name"
done

if [ "$mode" = "test" ] && [ -z "$selected_install_conf" ]; then
    selected_install_conf="conf/docker_test_settings.txt"
fi
install_conf="$selected_install_conf"

if [ "$mode" = "production" ] && [ -z "$install_conf" ]; then
    echo "Production deployments require --install_conf or --install_conf_map app,<path>."
    exit 1
fi

if [ -z "$skip_demo_data" ]; then
    if [ "$mode" = "test" ]; then
        skip_demo_data=false
    else
        skip_demo_data=true
    fi
fi

if [ -z "$skip_test_data" ]; then
    if [ "$mode" = "test" ]; then
        skip_test_data=false
    else
        skip_test_data=true
    fi
fi

if [ "$action" = "upgrade" ]; then
    skip_demo_data=true
    skip_test_data=true
fi

if [ ! -f "$compose_file" ]; then
    echo "Compose file '$compose_file' not found"
    exit 1
fi

if [ ! -f "$install_conf" ]; then
    echo "Install configuration '$install_conf' not found"
    exit 1
fi

repo_root="$(pwd)"
build_context_dir="$repo_root"
if [ ! -d "$build_context_dir" ]; then
    echo "Build context directory '$build_context_dir' not found"
    exit 1
fi

temp_install_conf=""
if [[ "$install_conf" = /* ]] && [[ "$install_conf" != "$build_context_dir/"* ]]; then
    temp_install_conf="$build_context_dir/.tmp_docker_install_conf_app_$$.txt"
    echo "Copying $install_conf into temporary file $temp_install_conf for Docker build/runtime."
    cp "$install_conf" "$temp_install_conf"
    install_conf="$temp_install_conf"
    cleanup_temp_conf() {
        if [ -n "$temp_install_conf" ] && [ -f "$temp_install_conf" ]; then
            rm -f "$temp_install_conf"
        fi
    }
    trap cleanup_temp_conf EXIT
fi

if [[ "$install_conf" = "$build_context_dir/"* ]]; then
    install_conf_container="${install_conf#$build_context_dir/}"
else
    install_conf_container="$install_conf"
fi

host_install_conf_path="$install_conf"
if [[ "$host_install_conf_path" != /* ]]; then
    host_install_conf_path="$repo_root/$host_install_conf_path"
fi

read_install_conf_value() {
    local key="$1"
    local file="$2"

    bash -c '
        set -a
        . "$1"
        key="$2"
        printf "%s" "${!key-}"
    ' _ "$file" "$key"
}

config_value_or_default() {
    local key="$1"
    local default_value="$2"
    local config_value=""
    local env_value="${!key:-}"

    if [ -n "$env_value" ]; then
        echo "$env_value"
        return 0
    fi

    config_value="$(read_install_conf_value "$key" "$host_install_conf_path")"
    if [ -n "$config_value" ]; then
        echo "$config_value"
    else
        echo "$default_value"
    fi
}

write_compose_env_file() {
    if [ "$mode" != "production" ]; then
        return 0
    fi

    cat > "$compose_env_file" << EOF
# Generated by container_install.sh from $install_conf_container.
# Used by Docker Compose/Podman Compose for docker-compose.prod.yml interpolation.
INSTALL_TYPE=dep
GIT_REVISION=$git_revision
INSTALL_CONF=$install_conf_container
INSTALL_PATH=$install_path
APACHE_CONF_PATH=$apache_conf_path
DJANGO_SETTINGS_PATH=$django_settings_path
APP_UID=$app_uid
APP_GID=$app_gid
APP_SHELL=$app_shell
APP_PORT=$app_port
DJANGO_DEBUG=$django_debug
DB_CONN_MAX_AGE=$db_conn_max_age
WEB_CONCURRENCY=$web_concurrency
GUNICORN_THREADS=$gunicorn_threads
GUNICORN_TIMEOUT=$gunicorn_timeout
GUNICORN_KEEPALIVE=$gunicorn_keepalive
SERVER_STATUS_SERVER_NAME=$apache_status_server_name
SERVER_STATUS_ALIASES=$apache_status_aliases
SERVER_STATUS_ALLOW_FROM=$apache_status_allow_from
APACHE_FORWARDED_PROTO=$apache_forwarded_proto
APACHE_FORWARDED_PORT=$apache_forwarded_port
EOF

    echo "Wrote Compose environment file: $compose_env_file"
}

compose_with_env_exec() {
    if [ "$mode" = "production" ] && [ -f "$compose_env_file" ]; then
        compose_exec --env-file "$compose_env_file" "$@"
    else
        compose_exec "$@"
    fi
}

set_engine

app_repo_path="${APP_REPO_PATH:-/srv/iskylims}"
config_install_path="$(read_install_conf_value "INSTALL_PATH" "$host_install_conf_path")"
install_path="${config_install_path:-/opt/iskylims}"
config_apache_conf_path="$(read_install_conf_value "APACHE_CONF_PATH" "$host_install_conf_path")"
apache_conf_path="${APACHE_CONF_PATH:-${config_apache_conf_path:-}}"
if [ -z "$apache_conf_path" ]; then
    apache_conf_path="$install_path/conf"
fi
config_django_settings_path="$(read_install_conf_value "DJANGO_SETTINGS_PATH" "$host_install_conf_path")"
django_settings_path="$(normalize_settings_bind_path "${DJANGO_SETTINGS_PATH:-${config_django_settings_path:-}}")"
app_uid="$(config_value_or_default APP_UID 1212)"
app_gid="$(config_value_or_default APP_GID 1212)"
app_shell="$(config_value_or_default APP_SHELL /sbin/nologin)"
app_port="$(config_value_or_default APP_PORT 8001)"
django_debug="$(config_value_or_default DJANGO_DEBUG false)"
db_conn_max_age="$(config_value_or_default DB_CONN_MAX_AGE 60)"
web_concurrency="$(config_value_or_default WEB_CONCURRENCY 2)"
gunicorn_threads="$(config_value_or_default GUNICORN_THREADS 2)"
gunicorn_timeout="$(config_value_or_default GUNICORN_TIMEOUT 300)"
gunicorn_keepalive="$(config_value_or_default GUNICORN_KEEPALIVE 5)"
config_dns_url="$(read_install_conf_value "DNS_URL" "$host_install_conf_path")"
apache_server_name="$(normalize_apache_server_name "${APACHE_SERVER_NAME:-${config_dns_url:-localhost}}")"
apache_log_name="$(printf '%s' "$apache_server_name" | tr -c 'A-Za-z0-9._-' '_' | sed 's/_$//')"
apache_status_server_name="$(normalize_apache_server_name "$(config_value_or_default SERVER_STATUS_SERVER_NAME "$apache_server_name")")"
apache_status_aliases="$(config_value_or_default SERVER_STATUS_ALIASES "127.0.0.1 localhost")"
apache_status_allow_from="$(config_value_or_default SERVER_STATUS_ALLOW_FROM "127.0.0.1 localhost")"
apache_forwarded_proto="$(config_value_or_default APACHE_FORWARDED_PROTO http)"
apache_forwarded_port="$(config_value_or_default APACHE_FORWARDED_PORT 8081)"
compose_env_file="$repo_root/.env.prod.file"
app_container=""
local_head_hash=""
local_head_short=""
app_image_name="${APP_IMAGE_NAME:-iskylims_app}"
image_id_before_build=""
image_id_after_build=""

# Check if a service exists in the compose file
#
# Parameters:
#   $1 - Service name to check
#
# Returns:
#   0 if the service exists, 1 otherwise
service_exists() {
    compose_with_env_exec -f "$compose_file" ps --services 2>/dev/null | grep -Fxq "$1"
}

# Return the name of the container for a given service name.
# The container name is different based on whether we are in test mode or not.
#
# Parameters:
#   $1 - Service name to return the container name for
#
# Returns:
#   The name of the container for the given service name
service_container_name() {
    local service_name="$1"
    if [ "$mode" = "test" ]; then
        case "$service_name" in
            db) echo "db" ;;
            samba) echo "samba" ;;
            app) echo "iskylims_app" ;;
            *) echo "" ;;
        esac
    else
        case "$service_name" in
            app) echo "iskylims_app" ;;
            samba) echo "samba" ;;
            *) echo "" ;;
        esac
    fi
}

# Resolve the container ID for the target app service.
#
# Returns:
#   Sets global variable `app_container` to a valid container name/ID.
#
# Errors:
#   Exits if unable to resolve a container for `app_service`.
resolve_app_container() {
    local container_name
    container_name="$(service_container_name "$app_service")"

    if [ -n "$container_name" ] && engine_exec inspect -f '{{.Id}}' "$container_name" >/dev/null 2>&1; then
        app_container="$container_name"
    else
        app_container="$(engine_exec ps -a --filter "label=com.docker.compose.service=${app_service}" --format '{{.ID}}' | head -n 1)"
    fi

    if [ -z "$app_container" ]; then
        echo "Error: unable to resolve container ID for service '$app_service'." >&2
        exit 1
    fi
}

try_resolve_app_container() {
    local container_name
    container_name="$(service_container_name "$app_service")"

    if [ -n "$container_name" ] && engine_exec inspect -f '{{.Id}}' "$container_name" >/dev/null 2>&1; then
        app_container="$container_name"
        return 0
    fi

    app_container="$(engine_exec ps -a --filter "label=com.docker.compose.service=${app_service}" --format '{{.ID}}' | head -n 1)"
    [ -n "$app_container" ]
}

# Ensure target app service container exists and is running.
#
# Errors:
#   Exits if container does not exist or is not running.
ensure_app_running() {
    resolve_app_container
    if ! engine_exec inspect -f '{{.State.Running}}' "$app_container" >/dev/null 2>&1; then
        echo "Error: service '$app_service' container does not exist."
        exit 1
    fi
    if [ "$(engine_exec inspect -f '{{.State.Running}}' "$app_container")" != "true" ]; then
        echo "Error: service '$app_service' container is not running. Showing logs:"
        engine_exec logs --tail 200 "$app_container"
        exit 1
    fi
}

print_local_source_diagnostics() {
    echo "Local source diagnostics:"
    if command -v git >/dev/null 2>&1 && git -C "$repo_root" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
        local_head_hash="$(git -C "$repo_root" rev-parse HEAD)"
        local_head_short="$(git -C "$repo_root" rev-parse --short HEAD)"
        echo "  local HEAD: $(git -C "$repo_root" log -1 --oneline)"
        echo "  local HEAD hash: $local_head_hash"
    else
        echo "  local git metadata unavailable"
    fi
}

print_existing_artifact_diagnostics() {
    echo "Image diagnostics before build:"
    if engine_exec image inspect "$app_image_name" >/dev/null 2>&1; then
        image_id_before_build="$(engine_exec image inspect -f '{{.Id}}' "$app_image_name" 2>/dev/null || true)"
        echo "  image before build: $image_id_before_build"
    else
        image_id_before_build=""
        echo "  image before build: none"
    fi
}

print_image_after_build() {
    echo "Image diagnostics after build:"
    if engine_exec image inspect "$app_image_name" >/dev/null 2>&1; then
        image_id_after_build="$(engine_exec image inspect -f '{{.Id}}' "$app_image_name" 2>/dev/null || true)"
        echo "  image after build: $image_id_after_build"
        if [ -n "$image_id_before_build" ] && [ "$image_id_before_build" = "$image_id_after_build" ]; then
            echo "  image id check: unchanged"
        elif [ -n "$image_id_before_build" ] && [ "$image_id_before_build" != "$image_id_after_build" ]; then
            echo "  image id check: changed"
        else
            echo "  image id check: created"
        fi
    else
        echo "  image after build: not found"
    fi
}

print_container_source_diagnostics() {
    local label="$1"
    local container_repo_head_hash=""
    local container_repo_head_short=""
    echo "$label"
    container_repo_head_hash="$(engine_exec exec "$app_container" sh -lc "
        if [ -d '$app_repo_path/.git' ]; then
            cd '$app_repo_path' && git rev-parse HEAD
        fi
    " 2>/dev/null | tail -n 1)"
    container_repo_head_short="$(engine_exec exec "$app_container" sh -lc "
        if [ -d '$app_repo_path/.git' ]; then
            cd '$app_repo_path' && git rev-parse --short HEAD
        fi
    " 2>/dev/null | tail -n 1)"
    if [ -n "$container_repo_head_hash" ]; then
        echo "  container /srv HEAD hash: $container_repo_head_hash"
    fi
    if [ -n "$local_head_hash" ] && [ -n "$container_repo_head_hash" ]; then
        if [ "$local_head_hash" = "$container_repo_head_hash" ]; then
            echo "  HEAD check: OK local=$local_head_short container=$container_repo_head_short"
        else
            echo "  HEAD check: MISMATCH local=$local_head_short container=$container_repo_head_short"
        fi
    fi
    engine_exec exec "$app_container" sh -lc "
        echo '  /srv/iskylims HEAD:'
        if [ -d '$app_repo_path/.git' ]; then
            cd '$app_repo_path' && git log -1 --oneline
        else
            echo 'not a git checkout'
        fi
    " || true
}

prepare_app_mount_permissions() {
    if [ "$mode" != "production" ]; then
        return 0
    fi

    echo "Preparing writable app mount permissions..."
    engine_exec exec --user 0 "$app_container" sh -lc "
        set -e
        mkdir -p '$install_path/logs' '$install_path/static' '$install_path/documents' '$install_path/cron' '$install_path/tmp'
        chown -R '$app_uid:$app_gid' '$install_path/logs' '$install_path/static' '$install_path/documents' '$install_path/cron' '$install_path/tmp'
        chmod -R u+rwX,g+rwX '$install_path/logs' '$install_path/static' '$install_path/documents'
        chmod 700 '$install_path/cron' '$install_path/tmp'
        if [ -f '$install_path/iskylims/settings.py' ]; then
            chown '$app_uid:$app_gid' '$install_path/iskylims/settings.py'
            chmod 0664 '$install_path/iskylims/settings.py'
        fi
        chmod -R o+rX '$install_path/static'
    "
}

prepare_host_bind_mount_permissions() {
    if [ "$mode" != "production" ]; then
        return 0
    fi

    local apache_conf_file
    local django_settings_dir

    echo "Preparing host bind mount permissions..."
    chmod_with_podman_fallback 0755 "$apache_conf_path"

    chown_with_podman_fallback "$app_uid:$app_gid" "/var/log/local/relecov-iskylims/apps"
    chmod_with_podman_fallback 0775 "/var/log/local/relecov-iskylims/apps"

    # UBI httpd runs as uid 1001 and group 0. This keeps the Apache log bind
    # writable without relying on Podman's :U ownership mutation.
    chown_with_podman_fallback "1001:0" "/var/log/local/iskylims/apache"
    chmod_with_podman_fallback 0775 "/var/log/local/iskylims/apache"

    if [ -f "$django_settings_path" ]; then
        django_settings_dir="$(dirname "$django_settings_path")"
        chmod_with_podman_fallback 0755 "$django_settings_dir"
        chown_with_podman_fallback "$app_uid:$app_gid" "$django_settings_path"
        chmod_with_podman_fallback 0664 "$django_settings_path"
    fi

    for apache_conf_file in \
        "$apache_conf_path/iskylims_apache_reverse_proxy.conf" \
        "$apache_conf_path/iskylims_apache_logs.conf" \
        "$apache_conf_path/iskylims_apache_server-status.conf"; do
        if [ -f "$apache_conf_file" ]; then
            chmod_with_podman_fallback 0664 "$apache_conf_file"
        fi
    done
}

# Remove stale test containers left over from previous runs.
#
# This function will only be executed in "test" mode when the engine is "podman".
# It only removes known test container names and never removes volumes.
cleanup_stale_test_containers() {
    if [ "$mode" != "test" ] || [ "$engine" != "podman" ]; then
        return 0
    fi

    local svc cname cstate
    for svc in db app samba; do
        cname="$(service_container_name "$svc")"
        if [ -z "$cname" ]; then
            continue
        fi
        if engine_exec inspect -f '{{.Id}}' "$cname" >/dev/null 2>&1; then
            cstate="$(engine_exec inspect -f '{{.State.Status}}' "$cname" 2>/dev/null || true)"
            if [ "$cstate" != "running" ]; then
                echo "Removing stale test container '$cname' (state: ${cstate:-unknown})"
                engine_exec rm -f "$cname" >/dev/null 2>&1 || true
            fi
        fi
    done
}

cleanup_stale_test_containers

if [ "$action" = "fix-permissions" ]; then
    echo "Repairing production container bind mount and volume permissions..."
    if ! mkdir -p "$apache_conf_path" "/var/log/local/iskylims/apache" "/var/log/local/relecov-iskylims/apps"; then
        echo "Error: unable to create required host bind/log directories. Check APACHE_CONF_PATH and log directory permissions." >&2
        exit 1
    fi
    prepare_django_settings_bind_mount "$django_settings_path"
    prepare_host_bind_mount_permissions
    write_compose_env_file
    if try_resolve_app_container && [ "$(engine_exec inspect -f '{{.State.Running}}' "$app_container" 2>/dev/null)" = "true" ]; then
        prepare_app_mount_permissions
        echo "Done repairing host bind mounts and mounted app volumes."
    else
        echo "Host bind mount permissions repaired."
        echo "The app container is not running, so named volumes were not repaired."
        echo "Start containers with Compose, then rerun this action to repair mounted app volumes."
    fi
    exit 0
fi

print_local_source_diagnostics
print_existing_artifact_diagnostics
echo "Deploying containers (compose file: $compose_file) with a pre-staged app image and GIT_REVISION=$git_revision..."
if ! mkdir -p "$apache_conf_path" "/var/log/local/iskylims/apache" "/var/log/local/relecov-iskylims/apps"; then
    echo "Error: unable to create required host bind/log directories. Check APACHE_CONF_PATH and log directory permissions." >&2
    exit 1
fi
prepare_django_settings_bind_mount "$django_settings_path"
if [ -f "$repo_root/conf/iskylims_apache_reverse_proxy.conf" ]; then
    render_apache_config \
        "$repo_root/conf/iskylims_apache_reverse_proxy.conf" \
        "$apache_conf_path/iskylims_apache_reverse_proxy.conf"
fi
if [ -f "$repo_root/conf/iskylims_apache_logs.conf" ]; then
    render_apache_config \
        "$repo_root/conf/iskylims_apache_logs.conf" \
        "$apache_conf_path/iskylims_apache_logs.conf"
fi
if [ -f "$repo_root/conf/iskylims_apache_server-status.conf" ]; then
    render_apache_config \
        "$repo_root/conf/iskylims_apache_server-status.conf" \
        "$apache_conf_path/iskylims_apache_server-status.conf"
fi
prepare_host_bind_mount_permissions
write_compose_env_file
INSTALL_TYPE="dep" GIT_REVISION="$git_revision" INSTALL_CONF="$install_conf_container" INSTALL_PATH="$install_path" APACHE_CONF_PATH="$apache_conf_path" DJANGO_SETTINGS_PATH="$django_settings_path" APP_UID="$app_uid" APP_GID="$app_gid" APP_SHELL="$app_shell" APP_PORT="$app_port" DJANGO_DEBUG="$django_debug" DB_CONN_MAX_AGE="$db_conn_max_age" WEB_CONCURRENCY="$web_concurrency" GUNICORN_THREADS="$gunicorn_threads" GUNICORN_TIMEOUT="$gunicorn_timeout" GUNICORN_KEEPALIVE="$gunicorn_keepalive" \
    compose_with_env_exec -f "$compose_file" build --no-cache \
    --build-arg INSTALL_TYPE="dep" \
    --build-arg GIT_REVISION="$git_revision" \
    --build-arg INSTALL_CONF="$install_conf_container" \
    --build-arg INSTALL_PATH="$install_path" \
    --build-arg APP_UID="$app_uid" \
    --build-arg APP_GID="$app_gid" \
    --build-arg APP_SHELL="$app_shell"
print_image_after_build
INSTALL_PATH="$install_path" APACHE_CONF_PATH="$apache_conf_path" DJANGO_SETTINGS_PATH="$django_settings_path" APP_UID="$app_uid" APP_GID="$app_gid" APP_SHELL="$app_shell" APP_PORT="$app_port" DJANGO_DEBUG="$django_debug" DB_CONN_MAX_AGE="$db_conn_max_age" WEB_CONCURRENCY="$web_concurrency" GUNICORN_THREADS="$gunicorn_threads" GUNICORN_TIMEOUT="$gunicorn_timeout" GUNICORN_KEEPALIVE="$gunicorn_keepalive" compose_with_env_exec -f "$compose_file" up -d

echo "Waiting 20 seconds for starting database and web services..."
sleep 20
ensure_app_running
print_container_source_diagnostics "Container diagnostics after startup:"
prepare_app_mount_permissions

container_install_conf_path="$install_conf_container"
if [[ "$container_install_conf_path" != /* ]]; then
    container_install_conf_path="$app_repo_path/$container_install_conf_path"
fi

if ! engine_exec exec -it "$app_container" test -f "$container_install_conf_path"; then
    echo "Copying install configuration into container at $container_install_conf_path"
    engine_exec cp "$host_install_conf_path" "${app_container}:$container_install_conf_path"
fi

script_args_before=""
if [ "$run_script_before" = true ]; then
    for val in "${migration_script_before[@]}"; do
        script_args_before+=" --script_before $(printf '%q' "$val")"
    done
fi

script_args_after=""
if [ "$run_script" = true ]; then
    for val in "${migration_script[@]}"; do
        script_args_after+=" --script_after $(printf '%q' "$val")"
    done
fi

if [ "$action" = "upgrade" ]; then
    echo "Running install.sh bootstrap inside the container (upgrade mode)"
    engine_exec exec -it "$app_container" bash -c "cd $app_repo_path && bash install.sh --bootstrap upgrade --git_revision \"$git_revision\" --conf \"$install_conf_container\" --tables --skip_apache_restart$script_args_before$script_args_after"
else
    echo "Running install.sh bootstrap inside the container (install mode)"
    engine_exec exec -it "$app_container" bash -c "cd $app_repo_path && bash install.sh --bootstrap install --git_revision \"$git_revision\" --conf \"$install_conf_container\" --skip_apache_restart$script_args_before$script_args_after"
fi

print_container_source_diagnostics "Container diagnostics after bootstrap:"

if ! engine_exec exec -it "$app_container" test -f "$install_path/manage.py"; then
    echo "Error: $install_path/manage.py not found after bootstrap. Showing logs:"
    engine_exec logs --tail 200 "$app_container"
    exit 1
fi

if [ "$skip_test_data" = false ]; then
    engine_exec exec -it "$app_container" python3 manage.py loaddata test/test_data.json
    engine_exec exec -it "$app_container" python3 manage.py shell -c "
from django.contrib.auth.models import Group, User
admin = User.objects.get(username='admin')
admin.groups.add(
    Group.objects.get(name='WetlabManager'),
    Group.objects.get(name='ServiceManager'),
)
print('admin groups:', list(admin.groups.values_list('name', flat=True)))
"
else
    echo "Skipping test data fixtures as requested"
fi

if [ "$skip_demo_data" = false ] && service_exists "samba"; then
    echo "Downloading and copying test files to the Samba container"
    if [ "$demo_data" == "false" ]; then
        wget https://zenodo.org/record/8091169/files/iskylims_demo_data.tar.gz
        demo_data="./iskylims_demo_data.tar.gz"
    fi
    engine_exec cp "$demo_data" samba:/mnt
    engine_exec exec -it samba tar -xf /mnt/iskylims_demo_data.tar.gz -C /mnt
    # Ensure extracted demo data can be traversed/read through SMB by non-owner users.
    engine_exec exec -it samba sh -lc '
        for root in /mnt/test_ngs_data /mnt/Runs; do
            if [ -d "$root" ]; then
                find "$root" -type d -exec chmod o+rx {} +
                find "$root" -type f -exec chmod o+r {} +
            fi
        done
    '

    echo "Deleting compressed test file"
    engine_exec exec -it samba rm /mnt/iskylims_demo_data.tar.gz

    if [ "$demo_data" == "false" ]; then
        rm -f "$demo_data"
    fi
else
    echo "Skipping Samba demo data load (flag enabled or service not present)"
fi

echo "Skipping crontab add/start (cron is managed by the container entrypoint)"

dns_url=""
local_ip=""
if [ -f "$host_install_conf_path" ]; then
    dns_url=$(grep -E "^DNS_URL=" "$host_install_conf_path" | tail -n 1 | cut -d= -f2- | sed "s/^['\"]//;s/['\"]$//")
    local_ip=$(grep -E "^LOCAL_SERVER_IP=" "$host_install_conf_path" | tail -n 1 | cut -d= -f2- | sed "s/^['\"]//;s/['\"]$//")
fi

access_urls=()
if [ -n "$dns_url" ] && [ "$dns_url" != "*" ]; then
    access_urls+=("http://${dns_url}:${app_port}")
fi
if [ -n "$local_ip" ] && [ "$local_ip" != "*" ]; then
    access_urls+=("http://${local_ip}:${app_port}")
fi
if [ ${#access_urls[@]} -eq 0 ]; then
    access_urls+=("http://localhost:${app_port}")
fi

echo "You can now access iSkyLIMS via: ${access_urls[*]}"
