#!/usr/bin/bash

container_install_script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "$container_install_script_dir/deployment/lib/container/common.sh"
# shellcheck disable=SC1091
source "$container_install_script_dir/deployment/lib/container/django.sh"

ISKYLIMS_VERSION="3.1.1"

# Print this application's stable outer-installer CLI and examples.
usage() {
cat << EOF
This script installs and upgrades the iskylims app.

Usage : $0 [--demo_data] [--git_revision] [--compose_file] [--install_conf] [--action] [--script] [--script_before] [--script_after] [--engine] [--test]
    Optional input data:
    --demo_data         | Provide already downloaded demo data from Zenodo
    --git_revision      | Specify the Git revision to install (default: main, or 'current' to use copied local sources)
    --compose_file      | Compose file to use (overrides default)
    --install_conf      | Single production settings file (host rendering, ephemeral build secret, and bootstrap)
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

# Supply iSkyLIMS Apache tokens. Token rendering and host-file installation are
# shared, but the set of virtual-host values is application-specific.
render_apache_config() {
    local src="$1"
    local dst="$2"
    render_config_template "$src" "$dst" 0664 \
        __ISKYLIMS_SERVER_NAME__ "$apache_server_name" \
        __ISKYLIMS_LOG_NAME__ "$apache_log_name" \
        __INSTALL_PATH__ "$install_path" \
        __APP_PORT__ "$app_port" \
        __ISKYLIMS_FORWARDED_PROTO__ "$apache_forwarded_proto" \
        __ISKYLIMS_FORWARDED_PORT__ "$apache_forwarded_port" \
        __GUNICORN_TIMEOUT__ "$gunicorn_timeout" \
        __SERVER_STATUS_SERVER_NAME__ "$apache_status_server_name" \
        __SERVER_STATUS_ALIASES__ "$apache_status_aliases" \
        __SERVER_STATUS_ALLOW_FROM__ "$apache_status_allow_from"
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

require_compose_file "$compose_file" || exit 1

repo_root="$(pwd)"
build_context_dir="$repo_root"
temp_install_conf=""
host_install_conf_path=""
prepare_install_configuration \
    "$install_conf" "$repo_root" "$build_context_dir" app "$mode" \
    host_install_conf_path install_conf_container temp_install_conf || exit 1
install_conf="$host_install_conf_path"
if [ -n "$temp_install_conf" ]; then
    # EXIT traps need a no-argument callback; deletion itself is shared.
    cleanup_temp_conf() { cleanup_files "$temp_install_conf"; }
    trap cleanup_temp_conf EXIT
fi

# Write this stack's complete Compose interpolation contract. The keys remain
# application-owned because they must match its production Compose file.
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
# Compose reads this host path to create an ephemeral build-secret mount. The
# value is used by Compose itself and is not stored in the resulting image.
INSTALL_CONF_HOST=$host_install_conf_path
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
django_settings_path="$(normalize_bind_file_path "${DJANGO_SETTINGS_PATH:-${config_django_settings_path:-}}" "$install_path/iskylims/settings.py" settings.py)"
app_uid="$(config_value_or_default APP_UID "$host_install_conf_path" 1212)"
app_gid="$(config_value_or_default APP_GID "$host_install_conf_path" 1212)"
app_shell="$(config_value_or_default APP_SHELL "$host_install_conf_path" /sbin/nologin)"
app_port="$(config_value_or_default APP_PORT "$host_install_conf_path" 8001)"
django_debug="$(config_value_or_default DJANGO_DEBUG "$host_install_conf_path" false)"
db_conn_max_age="$(config_value_or_default DB_CONN_MAX_AGE "$host_install_conf_path" 60)"
web_concurrency="$(config_value_or_default WEB_CONCURRENCY "$host_install_conf_path" 2)"
gunicorn_threads="$(config_value_or_default GUNICORN_THREADS "$host_install_conf_path" 2)"
gunicorn_timeout="$(config_value_or_default GUNICORN_TIMEOUT "$host_install_conf_path" 300)"
gunicorn_keepalive="$(config_value_or_default GUNICORN_KEEPALIVE "$host_install_conf_path" 5)"
config_dns_url="$(read_install_conf_value "DNS_URL" "$host_install_conf_path")"
apache_server_name="$(normalize_apache_server_name "${APACHE_SERVER_NAME:-${config_dns_url:-localhost}}")"
apache_log_name="$(printf '%s' "$apache_server_name" | tr -c 'A-Za-z0-9._-' '_' | sed 's/_$//')"
apache_status_server_name="$(normalize_apache_server_name "$(config_value_or_default SERVER_STATUS_SERVER_NAME "$host_install_conf_path" "$apache_server_name")")"
apache_status_aliases="$(config_value_or_default SERVER_STATUS_ALIASES "$host_install_conf_path" "127.0.0.1 localhost")"
apache_status_allow_from="$(config_value_or_default SERVER_STATUS_ALLOW_FROM "$host_install_conf_path" "127.0.0.1 localhost")"
apache_forwarded_proto="$(config_value_or_default APACHE_FORWARDED_PROTO "$host_install_conf_path" https)"
apache_forwarded_port="$(config_value_or_default APACHE_FORWARDED_PORT "$host_install_conf_path" 443)"
compose_env_file="$repo_root/.env.prod.file"
app_container=""
if [ "$mode" = "production" ]; then
    app_image_name="${APP_IMAGE_NAME:-relecov-iskylims:local}"
else
    app_image_name="${APP_IMAGE_NAME:-iskylims_app}"
fi

# Map services to this stack's explicit legacy container names.
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

# Repair iSkyLIMS-specific writable directories inside its mounted app volume.
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
        chmod -R o+rX '$install_path/static'
    "
    prepare_django_container_settings_permissions \
        "$app_container" "$install_path/iskylims/settings.py" "$app_uid" "$app_gid"
}

# Apply ownership required by iSkyLIMS and the UBI Apache runtime on host binds.
prepare_host_bind_mount_permissions() {
    if [ "$mode" != "production" ]; then
        return 0
    fi

    local django_settings_dir
    local -a host_permission_spec

    echo "Preparing host bind mount permissions..."
    host_permission_spec=(
        "$apache_conf_path|-|0755"
        "/var/log/local/relecov-iskylims/apps|$app_uid:$app_gid|0775"
        # UBI httpd runs as uid 1001 and group 0. This keeps the Apache log bind
        # writable without relying on Podman's :U ownership mutation.
        "/var/log/local/relecov-iskylims/apache|1001:0|0775"
        "$apache_conf_path/iskylims_apache_reverse_proxy.conf|-|0664"
        "$apache_conf_path/iskylims_apache_logs.conf|-|0664"
        "$apache_conf_path/iskylims_apache_server-status.conf|-|0664"
    )

    if [ -f "$django_settings_path" ]; then
        django_settings_dir="$(dirname "$django_settings_path")"
        host_permission_spec+=(
            "$django_settings_dir|-|0755"
            "$django_settings_path|$app_uid:$app_gid|0664"
        )
    fi

    apply_host_permission_spec "${host_permission_spec[@]}"
}

if [ "$action" = "fix-permissions" ]; then
    echo "Repairing production container bind mount and volume permissions..."
    if ! mkdir -p "$apache_conf_path" "/var/log/local/relecov-iskylims/apache" "/var/log/local/relecov-iskylims/apps"; then
        echo "Error: unable to create required host bind/log directories. Check APACHE_CONF_PATH and log directory permissions." >&2
        exit 1
    fi
    prepare_django_settings_bind_mount "$repo_root/conf/template_settings.txt" "$django_settings_path" "$host_install_conf_path"
    prepare_host_bind_mount_permissions
    write_compose_env_file
    if app_container="$(resolve_service_container "$app_service" 2>/dev/null)" \
        && [ "$(engine_exec inspect -f '{{.State.Running}}' "$app_container" 2>/dev/null)" = "true" ]; then
        prepare_app_mount_permissions
        echo "Done repairing host bind mounts and mounted app volumes."
    else
        echo "Host bind mount permissions repaired."
        echo "The app container is not running, so named volumes were not repaired."
        echo "Start containers with Compose, then rerun this action to repair mounted app volumes."
    fi
    exit 0
fi

expected_head_hash="$(repository_revision "$repo_root" full)" || true
expected_head_short="$(repository_revision "$repo_root" short)" || true
image_id_before_build="$(compose_service_image_id "$compose_file" "$app_service")"
print_prebuild_diagnostics \
    "Pre-build diagnostics" "$repo_root" "$git_revision" "$image_id_before_build"
echo "Deploying containers (compose file: $compose_file) with a pre-staged app image and GIT_REVISION=$git_revision..."
if ! mkdir -p "$apache_conf_path" "/var/log/local/relecov-iskylims/apache" "/var/log/local/relecov-iskylims/apps"; then
    echo "Error: unable to create required host bind/log directories. Check APACHE_CONF_PATH and log directory permissions." >&2
    exit 1
fi
prepare_django_settings_bind_mount "$repo_root/conf/template_settings.txt" "$django_settings_path" "$host_install_conf_path"
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
validate_compose_configuration "$compose_file" || exit 1
if [ "$mode" = "production" ]; then
    # Use the engine directly because older Compose schemas cannot express
    # build secrets. The explicit tag matches docker-compose.prod.yml.
    engine_build --no-cache \
        --secret "id=install_conf,src=$host_install_conf_path" \
        --build-arg INSTALL_TYPE="dep" \
        --build-arg GIT_REVISION="$git_revision" \
        --build-arg INSTALL_CONF="$install_conf_container" \
        --build-arg USE_INSTALL_CONF_SECRET="true" \
        --build-arg RENDER_DJANGO_SETTINGS="false" \
        --build-arg INSTALL_PATH="$install_path" \
        --build-arg APP_UID="$app_uid" \
        --build-arg APP_GID="$app_gid" \
        --build-arg APP_SHELL="$app_shell" \
        --tag "$app_image_name" \
        "$build_context_dir"
else
    INSTALL_TYPE="dep" GIT_REVISION="$git_revision" INSTALL_CONF="$install_conf_container" INSTALL_PATH="$install_path" \
        compose_with_env_exec -f "$compose_file" build --no-cache
fi
image_id_after_build="$(compose_service_image_id "$compose_file" "$app_service")"
print_image_after_diagnostics \
    "Image diagnostics after build:" "$image_id_before_build" "$image_id_after_build"
INSTALL_PATH="$install_path" APACHE_CONF_PATH="$apache_conf_path" DJANGO_SETTINGS_PATH="$django_settings_path" APP_UID="$app_uid" APP_GID="$app_gid" APP_SHELL="$app_shell" APP_PORT="$app_port" DJANGO_DEBUG="$django_debug" DB_CONN_MAX_AGE="$db_conn_max_age" WEB_CONCURRENCY="$web_concurrency" GUNICORN_THREADS="$gunicorn_threads" GUNICORN_TIMEOUT="$gunicorn_timeout" GUNICORN_KEEPALIVE="$gunicorn_keepalive" compose_with_env_exec -f "$compose_file" up -d

echo "Waiting 20 seconds for starting database and web services..."
sleep 20
app_container="$(resolve_service_container "$app_service")" || exit 1
ensure_service_running "$app_service" "$app_container" >/dev/null
print_container_repository_diagnostics \
    "Container diagnostics after startup:" "$app_container" "$app_repo_path" \
    "$expected_head_hash" "$expected_head_short"
prepare_app_mount_permissions

container_install_conf_path="$install_conf_container"
if [[ "$container_install_conf_path" != /* ]]; then
    container_install_conf_path="$app_repo_path/$container_install_conf_path"
fi

echo "Staging protected runtime configuration at $container_install_conf_path"
stage_container_runtime_config \
    "$app_container" "$host_install_conf_path" "$container_install_conf_path" \
    "$app_uid" "$app_gid" || exit 1

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

# Bootstrap consumes a temporary runtime copy. Remove that exact file so the
# operator configuration remains only on the protected host after deployment.
if [ "$mode" = "production" ]; then
    remove_container_runtime_config \
        "$app_container" "$container_install_conf_path" || exit 1
fi

print_container_repository_diagnostics \
    "Container diagnostics after bootstrap:" "$app_container" "$app_repo_path" \
    "$expected_head_hash" "$expected_head_short"

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

run_standard_smoke_test \
    "$repo_root/scripts/smoke_test.sh" "$mode" "$engine" "$compose_file" \
    "${compose_env_file:-}" || exit 1

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
