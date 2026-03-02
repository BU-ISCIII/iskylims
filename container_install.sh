#!/usr/bin/bash

ISKYLIMS_VERSION="3.1.0"

usage() {
cat << EOF
This script installs and upgrades the iskylims app.

Usage : $0 [--demo_data] [--install_type] [--git_revision] [--compose_file] [--install_conf] [--action] [--script] [--script_before] [--script_after] [--engine] [--test]
    Optional input data:
    --demo_data         | Provide already downloaded demo data from Zenodo
    --install_type      | Specify the installation type for iSkyLIMS (default: full)
    --git_revision      | Specify the Git revision to install (default: main)
    --compose_file      | Compose file to use (overrides default)
    --install_conf      | Settings file consumed during container image build (mandatory for production)
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

    Upgrade an existing production deployment using the same database:
    bash $0 --install_conf conf/my_prod_settings.txt --action upgrade

    Install demo container system with local services
    bash $0 --test

    Provide already downloaded data from Zenodo (compressed) for test environment
    bash $0 --demo_data /path/to/iskylims_demo_data.tar.gz

    Speficy a custom installation using the Git revision "develop":
    bash $0 --install_type app --git_revision develop

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
        --install_type)      set -- "$@" -i ;;
        --git_revision)      set -- "$@" -g ;;
        --compose_file)      set -- "$@" -c ;;
        --install_conf)      set -- "$@" -s ;;
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
options=":d:i:g:c:s:a:m:b:f:e:vhntp"
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
    if [ -z "$install_conf" ]; then
        install_conf="conf/docker_test_settings.txt"
    fi
else
    if [ -z "$compose_file" ]; then
        compose_file="docker-compose.prod.yml"
    fi
fi

if [ "$mode" = "production" ] && [ -z "$install_conf" ]; then
    echo "Production deployments require --install_conf pointing to your settings file."
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
temp_install_conf=""
if [[ "$install_conf" = /* ]] && [[ "$install_conf" != "$repo_root/"* ]]; then
    temp_install_conf="$repo_root/.tmp_docker_install_conf_$$.txt"
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

if [[ "$install_conf" = "$repo_root/"* ]]; then
    install_conf_container="${install_conf#$repo_root/}"
else
    install_conf_container="$install_conf"
fi

set_engine

app_service="${APP_SERVICE:-app}"
app_repo_path="${APP_REPO_PATH:-/srv/iskylims}"
app_install_path="${APP_INSTALL_PATH:-/opt/iskylims}"
app_port="${APP_PORT:-8001}"
app_container=""

service_exists() {
    compose_exec -f "$compose_file" ps --services 2>/dev/null | grep -Fxq "$1"
}

resolve_app_container() {
    app_container="$(compose_exec -f "$compose_file" ps -q "$app_service" | head -n 1)"
    if [ -z "$app_container" ]; then
        echo "Error: unable to resolve container ID for service '$app_service'."
        exit 1
    fi
}

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

echo "Deploying containers (compose file: $compose_file) with INSTALL_TYPE=dep and GIT_REVISION=$git_revision..."
INSTALL_TYPE="dep" GIT_REVISION="$git_revision" INSTALL_CONF="$install_conf_container" \
    compose_exec -f "$compose_file" build --no-cache \
    --build-arg INSTALL_TYPE="dep" \
    --build-arg GIT_REVISION="$git_revision" \
    --build-arg INSTALL_CONF="$install_conf_container"
compose_exec -f "$compose_file" up -d

echo "Waiting 20 seconds for starting database and web services..."
sleep 20
ensure_app_running

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

if ! engine_exec exec -it "$app_container" test -f "$app_install_path/manage.py"; then
    echo "Error: $app_install_path/manage.py not found after install.sh. Showing logs:"
    engine_exec logs --tail 200 "$app_container"
    exit 1
fi

if [ "$skip_test_data" = false ]; then
    engine_exec exec -it "$app_container" python3 manage.py loaddata test/test_data.json
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
