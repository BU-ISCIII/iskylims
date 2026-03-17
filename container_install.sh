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
    --action            | install (default) or upgrade, to control DB initialisation steps
    --script            | Run a Django migration script after migrations (can be repeated)
    --script_before     | Run a Django migration script before migrations (can be repeated)
    --script_after      | Run a Django migration script after migrations (can be repeated)
    --skip_demo_data    | Skip downloading/copying demo data to samba container
    --skip_test_data    | Skip loading test fixtures (test/test_data.json)
    --engine            | Container engine to use: docker (default) or podman
    --test              | Use development/test compose file and sample data

Examples:
    Deploy production container pointing to an external DB/Samba:
    bash $0 --install_conf conf/my_prod_settings.txt

    Deploy production with service-specific settings mapping:
    bash $0 --install_conf_map app,conf/docker_production_settings.txt

    Upgrade an existing production deployment using the same database:
    bash $0 --install_conf conf/my_prod_settings.txt --action upgrade

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
            if [[ "$action" != "install" && "$action" != "upgrade" ]]; then
                echo "Invalid action '$action'. Use install or upgrade."
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

set_engine

app_repo_path="${APP_REPO_PATH:-/srv/iskylims}"
app_install_path="${APP_INSTALL_PATH:-/opt/iskylims}"
app_port="${APP_PORT:-8001}"
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
    compose_exec -f "$compose_file" ps --services 2>/dev/null | grep -Fxq "$1"
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

print_local_source_diagnostics
print_existing_artifact_diagnostics
echo "Deploying containers (compose file: $compose_file) with INSTALL_TYPE=dep and GIT_REVISION=$git_revision..."
INSTALL_TYPE="dep" GIT_REVISION="$git_revision" INSTALL_CONF="$install_conf_container" \
    compose_exec -f "$compose_file" build --no-cache \
    --build-arg INSTALL_TYPE="dep" \
    --build-arg GIT_REVISION="$git_revision" \
    --build-arg INSTALL_CONF="$install_conf_container"
print_image_after_build
compose_exec -f "$compose_file" up -d

echo "Waiting 20 seconds for starting database and web services..."
sleep 20
ensure_app_running
print_container_source_diagnostics "Container diagnostics after startup:"

app_uid="${APP_UID:-1212}"
app_gid="${APP_GID:-1212}"
echo "Ensuring runtime directories are writable by ${app_uid}:${app_gid}"
engine_exec exec -u 0 -it "$app_container" sh -lc "mkdir -p ${app_install_path}/documents ${app_install_path}/logs ${app_install_path}/static ${app_install_path}/cron ${app_install_path}/tmp && chown -R ${app_uid}:${app_gid} ${app_install_path}/documents ${app_install_path}/logs ${app_install_path}/static ${app_install_path}/cron ${app_install_path}/tmp"

host_install_conf_path="$install_conf"
if [[ "$host_install_conf_path" != /* ]]; then
    host_install_conf_path="$repo_root/$host_install_conf_path"
fi

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
    echo "Running install.sh upgrade inside the container"
    engine_exec exec -it "$app_container" bash -c "cd $app_repo_path && bash install.sh --upgrade app --git_revision \"$git_revision\" --conf \"$install_conf_container\" --skip_apache_restart$script_args_before$script_args_after"
else
    echo "Running install.sh install inside the container"
    engine_exec exec -it "$app_container" bash -c "cd $app_repo_path && bash install.sh --install app --git_revision \"$git_revision\" --conf \"$install_conf_container\" --skip_apache_restart$script_args_before$script_args_after"
fi

print_container_source_diagnostics "Container diagnostics after install.sh:"

if ! engine_exec exec -it "$app_container" test -f "$app_install_path/manage.py"; then
    echo "Error: $app_install_path/manage.py not found after install.sh. Showing logs:"
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
