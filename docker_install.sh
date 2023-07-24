#!/usr/bin/bash

ISKYLIMS_VERSION="3.0.0"

usage() {
cat << EOF
This script install and upgrade the iskylims app.

usage : $0 --demo_data
    Optional input data:
    --demo_data  | provide already dowloaded demo data from zenodo


Examples:
    Install demo docker system
    bash $0 

    Provide already downloaded data from zenodo (compressed)
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
        --demo_data)     set -- "$@" -d ;;

    # ADITIONAL
        --help)          set -- "$@" -h ;;
        --version)       set -- "$@" -v ;;
    # PASSING VALUE IN PARAMETER
        *)               set -- "$@" "$arg" ;;
    esac
done

# SETTING DEFAULT VALUES
demo_data=false

# PARSE VARIABLE ARGUMENTS WITH getops
options=":d:vh"
while getopts $options opt; do
    case $opt in
        d)
            demo_data=$OPTARG
            ;;
        h )
            usage
            exit 1
            ;;
        v )
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

echo "Deploying test containers..."
docker compose build --no-cache 
docker compose up -d

echo "Waiting 20 seconds for starting database and web services..."
sleep 20

echo "Creating the database structure for iSkyLIMS"
docker exec -it iskylims_app python3 manage.py migrate
docker exec -it iskylims_app python3 manage.py makemigrations django_utils core wetlab drylab
docker exec -it iskylims_app python3 manage.py migrate

echo "Loading in database initial data"
docker exec -it iskylims_app python3 manage.py loaddata conf/first_install_tables.json

echo "Creating super user "
docker exec -it iskylims_app python3 manage.py createsuperuser

echo "Download testing files and copy it to samba container"
if [ "$demo_data" == "false" ];then
    wget https://zenodo.org/record/8091169/files/iskylims_demo_data.tar.gz
    demo_data="./iskylims_demo_data.tar.gz"
fi
docker cp $demo_data samba:/mnt
docker exec -it samba tar -xf  /mnt/iskylims_demo_data.tar.gz -C /mnt

echo "deleting compress testing file"
docker exec -it samba rm /mnt/iskylims_demo_data.tar.gz

if [ "$demo_data" == "false" ];then
    rm -f $demo_data
fi

echo "Running crontab"
docker exec -it iskylims_app python3 manage.py crontab add
echo "You can now access iskylims via: http://localhost:8001"