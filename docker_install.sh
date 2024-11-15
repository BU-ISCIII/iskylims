#!/usr/bin/bash

ISKYLIMS_VERSION="3.0.0"

usage() {
cat << EOF
This script installs and upgrades the iskylims app.

Usage : $0 [--demo_data] [--install_type] [--git_revision]
    Optional input data:
    --demo_data         | Provide already downloaded demo data from Zenodo
    --install_type      | Specify the installation type for iSkyLIMS (default: full)
    --git_revision      | Specify the Git revision to install (default: main)

Examples:
    Install demo docker system
    bash $0 

    Provide already downloaded data from Zenodo (compressed)
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

        # ADDITIONAL
        --help)              set -- "$@" -h ;;
        --version)           set -- "$@" -v ;;
        # PASSING VALUE IN PARAMETER
        *)                   set -- "$@" "$arg" ;;
    esac
done

# SETTING DEFAULT VALUES
demo_data=false
install_type="app"
git_revision="main"

# PARSE VARIABLE ARGUMENTS WITH getopts
options=":d:i:g:vh"
while getopts $options opt; do
    case $opt in
        d)
            demo_data=$OPTARG
            ;;
        i)
            install_type=$OPTARG
            ;;
        g)
            git_revision=$OPTARG
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

echo "Deploying test containers with INSTALL_TYPE=$install_type and GIT_REVISION=$git_revision..."
docker compose build --no-cache --build-arg INSTALL_TYPE=$install_type --build-arg GIT_REVISION=$git_revision
docker compose up -d

echo "Waiting 20 seconds for starting database and web services..."
sleep 20

echo "Creating the database structure for iSkyLIMS"
docker exec -it iskylims_app python3 manage.py migrate
docker exec -it iskylims_app python3 manage.py makemigrations django_utils core wetlab drylab
docker exec -it iskylims_app python3 manage.py migrate

echo "Creating super user"
docker exec -it iskylims_app python3 manage.py createsuperuser

echo "Loading initial data into the database"
docker exec -it iskylims_app python3 manage.py loaddata conf/first_install_tables.json
docker exec -it iskylims_app python3 manage.py loaddata test/test_data.json

echo "Downloading and copying test files to the Samba container"
if [ "$demo_data" == "false" ]; then
    wget https://zenodo.org/record/8091169/files/iskylims_demo_data.tar.gz
    demo_data="./iskylims_demo_data.tar.gz"
fi
docker cp $demo_data samba:/mnt
docker exec -it samba tar -xf /mnt/iskylims_demo_data.tar.gz -C /mnt

echo "Deleting compressed test file"
docker exec -it samba rm /mnt/iskylims_demo_data.tar.gz

if [ "$demo_data" == "false" ]; then
    rm -f $demo_data
fi

echo "Running crontab"
docker exec -it iskylims_app python3 manage.py crontab add
docker exec -it iskylims_app service cron start

echo "You can now access iSkyLIMS via: http://localhost:8001"
