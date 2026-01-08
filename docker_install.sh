#!/usr/bin/bash

ISKYLIMS_VERSION="3.0.0"

usage() {
cat << EOF
This script installs and upgrades the iskylims app.

Usage : $0 [--demo_data] [--install_type] [--git_revision] [--compose_file] [--install_conf] [--action] [--script] [--test]
    Optional input data:
    --demo_data         | Provide already downloaded demo data from Zenodo
    --install_type      | Specify the installation type for iSkyLIMS (default: full)
    --git_revision      | Specify the Git revision to install (default: main)
    --compose_file      | docker compose file to use (overrides default)
    --install_conf      | Settings file consumed during docker image build (mandatory for production)
    --action            | install (default) or upgrade, to control DB initialisation steps
    --script            | Run a Django migration script (can be repeated)
    --skip_demo_data    | Skip downloading/copying demo data to samba container
    --skip_test_data    | Skip loading test fixtures (test/test_data.json)
    --test              | Use development/test compose file and sample data

Examples:
    Deploy production container pointing to an external DB/Samba:
    bash $0 --install_conf conf/my_prod_settings.txt

    Upgrade an existing production deployment using the same database:
    bash $0 --install_conf conf/my_prod_settings.txt --action upgrade

    Install demo docker system with local services
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
        --skip_demo_data)    set -- "$@" -n ;;
        --skip_test_data)    set -- "$@" -t ;;
        --test)              set -- "$@" -p ;;

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
skip_demo_data=""
skip_test_data=""
mode="production"
action="install"
run_script=false
migration_script=()

# PARSE VARIABLE ARGUMENTS WITH getopts
options=":d:i:g:c:s:a:m:vhntp"
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
        compose_file="docker-compose.yml"
    fi
    if [ -z "$install_conf" ]; then
        install_conf="conf/docker_install_settings.txt"
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
if [[ "$install_conf" = /* ]] && [[ "$install_conf" != "$repo_root/"* ]]; then
    install_conf_copy="conf/docker_runtime_settings.txt"
    echo "Copying $install_conf into repository as $install_conf_copy for Docker build/runtime."
    cp "$install_conf" "$install_conf_copy"
    install_conf="$install_conf_copy"
fi

service_exists() {
    docker compose -f "$compose_file" ps --services 2>/dev/null | grep -Fxq "$1"
}

ensure_app_running() {
    if ! docker inspect -f '{{.State.Running}}' iskylims_app >/dev/null 2>&1; then
        echo "Error: iskylims_app container does not exist."
        exit 1
    fi
    if [ "$(docker inspect -f '{{.State.Running}}' iskylims_app)" != "true" ]; then
        echo "Error: iskylims_app container is not running. Showing logs:"
        docker logs --tail 200 iskylims_app
        exit 1
    fi
}

echo "Deploying containers (compose file: $compose_file) with INSTALL_TYPE="dep" and GIT_REVISION=$git_revision..."
docker compose -f "$compose_file" build --no-cache --build-arg INSTALL_TYPE="dep" --build-arg GIT_REVISION="$git_revision" --build-arg INSTALL_CONF="$install_conf"
docker compose -f "$compose_file" up -d

echo "Waiting 20 seconds for starting database and web services..."
sleep 20
ensure_app_running

script_args=""
if [ "$run_script" = true ]; then
    for val in "${migration_script[@]}"; do
        script_args+=" --script $(printf '%q' "$val")"
    done
fi

if [ "$action" = "upgrade" ]; then
    echo "Running install.sh upgrade inside the container"
    docker exec -it iskylims_app bash -c "cd /srv/iskylims && bash install.sh --upgrade app --git_revision \"$git_revision\" --conf \"$install_conf\" --skip_apache_restart$script_args"
else
    echo "Running install.sh install inside the container"
    docker exec -it iskylims_app bash -c "cd /srv/iskylims && bash install.sh --install app --git_revision \"$git_revision\" --conf \"$install_conf\" --skip_apache_restart$script_args"
fi

if ! docker exec -it iskylims_app test -f /opt/iskylims/manage.py; then
    echo "Error: /opt/iskylims/manage.py not found after install.sh. Showing logs:"
    docker logs --tail 200 iskylims_app
    exit 1
fi

if [ "$skip_test_data" = false ]; then
    docker exec -it iskylims_app python3 manage.py loaddata test/test_data.json
else
    echo "Skipping test data fixtures as requested"
fi

if [ "$skip_demo_data" = false ] && service_exists "samba"; then
    echo "Downloading and copying test files to the Samba container"
    if [ "$demo_data" == "false" ]; then
        wget https://zenodo.org/record/8091169/files/iskylims_demo_data.tar.gz
        demo_data="./iskylims_demo_data.tar.gz"
    fi
    docker cp "$demo_data" samba:/mnt
    docker exec -it samba tar -xf /mnt/iskylims_demo_data.tar.gz -C /mnt

    echo "Deleting compressed test file"
    docker exec -it samba rm /mnt/iskylims_demo_data.tar.gz

    if [ "$demo_data" == "false" ]; then
        rm -f "$demo_data"
    fi
else
    echo "Skipping Samba demo data load (flag enabled or service not present)"
fi

echo "Running crontab"
docker exec -it iskylims_app python3 manage.py crontab add
docker exec -it iskylims_app service cron start

echo "You can now access iSkyLIMS via: http://localhost:8001"
