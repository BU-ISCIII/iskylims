#!/bin/bash

#=============================================================
# HEADER
#=============================================================

#INSTITUTION:ISCIII
#CENTRE:BU-ISCIII
#AUTHOR: Luis Chapado
SCRIPT_VERSION=0.1
#CREATED: 01 September 2021
#
#
#DESCRIPTION: 
#
#
#================================================================
# END_OF_HEADER
#================================================================

ISKYLIMS_VERSION=2.0.0

DB_USER=django
DB_PASS=django77
DB_SERVER_IP=localhost
DB_PORT=3306

db_check(){
    mysqladmin -h $DB_SERVER_IP -u$DB_USER -p$DB_PASS -P$DB_PORT processlist >null ###user should have mysql permission on remote server.

    if ! [ $? -eq 0 ]; then
        echo -e "${RED}ERROR : Unable to connect to database. Check if your database is running and accessible${NC}"
        exit 1
    fi
    RESULT=`mysqlshow --user=$DB_USER --password=$DB_PASS --host=$DB_SERVER_IP --port=$DB_PORT | grep -o iSkyLIMS`

    if  ! [ "$RESULT" == "iSkyLIMS" ] ; then
        echo -e "${RED}ERROR : iSkyLIMS database is not defined yet ${NC}"
        echo -e "${RED}ERROR : Create iSkyLIMS database on your mysql server and run again the installation script ${NC}"
        exit 1    
    fi
}

#================================================================
#SET TEMINAL COLORS
#================================================================
YELLOW='\033[0;33m'
WHITE='\033[0;37m'
CYAN='\033[0;36m'
BLUE='\033[0;34m'
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

#================================================================
# MAIN_BODY
#================================================================

printf "\n\n%s"
printf "${YELLOW}------------------${NC}\n"
printf "%s"
printf "${YELLOW}Starting iSkyLIMS Installation version: ${ISKYLIMS_VERSION}${NC}\n"
printf "%s"
printf "${YELLOW}------------------${NC}\n\n"

#================================================================
#CHECK REQUIREMENTS BEFORE STARTING INSTALLATION
#================================================================
# Check that script is run as root
if [[ $EUID -ne 0 ]]; then
    printf "\n\n%s"
    printf "${RED}------------------${NC}\n"
    printf "%s"
    printf "${RED}Exiting installation. This script must be run as root ${NC}\n" 
    printf "\n\n%s"
    printf "${RED}------------------${NC}\n"
    printf "%s"
    exit 1
fi

echo "Checking main requirements"
db_check
echo "Successful check for database"
linux_distribution=$(lsb_release -i | cut -f 2-)
#================================================================

read -p "Are you sure you want to install iSkyLIMS in this server? " -n 1 -r
echo    # (optional) move to a new line
if [[ ! $REPLY =~ ^[Yy]$ ]] ; then
    exit 1
fi

#================================================================
if [[ $linux_distribution == "Ubuntu" ]]; then 
    echo "Software installation for Ubuntu"
    apt-get update && apt-get upgrade -y
    apt-get install -y \
        lightdm git apt-utils libcairo2 libcairo2-dev  wget gnuplot python3-pip \
        libmysqlclient-dev apache2 apache2-dev vim
    #apt-get install build-essential  -y 
    #apt-get install libghc-zlib-dev libbz2-dev libssl1.0-dev -y
    #apt-get install git libpango1.0 libpango1.0-dev   -y
fi

#mkdir -p /opt/interop
#cd /opt/interop

#wget https://github.com/Illumina/interop/releases/download/v1.1.15/InterOp-1.1.15-Linux-GNU.tar.gz
#tar -xf  InterOp-1.1.15-Linux-GNU.tar.gz
#ln -s InterOp-1.1.15-Linux-GNU interop
#rm InterOp-1.1.15-Linux-GNU.tar.gz


echo "iSkyLIMS installation"
mkdir -p /opt/iSkyLIMS
cd /opt/iSkyLIMS
git clone https://github.com/BU-ISCIII/iSkyLIMS.git .

git submodule init
git checkout develop
git submodule init
git submodule update --checkout


mkdir -p /opt/iSkyLIMS/documents/wetlab/tmp
mkdir -p /opt/iSkyLIMS/documents/drylab
mkdir -p /opt/iSkyLIMS/logs

# install virtual environment
if [[ $linux_distribution == "Ubuntu" ]]; then 
    echo "Creating virtual environment"
    pip3 install virtualenv 
    virtualenv --python=/usr/bin/python3 virtualenv
    
    
fi

source virtualenv/bin/activate

# Starting iSkyLIMS
python3 -m pip install -r conf/pythonPackagesRequired.txt
django-admin startproject iSkyLIMS .
grep ^SECRET iSkyLIMS/settings.py > ~/.secret


# Copying config files and script
cp conf/docker_settings.py /opt/iSkyLIMS/iSkyLIMS/settings.py
cp conf/urls.py /opt/iSkyLIMS/iSkyLIMS/

sed -i "/^SECRET/c\\$(cat ~/.secret)" iSkyLIMS/settings.py

