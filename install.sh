#!/bin/bash

#=============================================================
# HEADER
#=============================================================

#INSTITUTION:ISCIII
#CENTRE:BU-ISCIII
#AUTHOR: Luis Chapado
SCRIPT_VERSION=0.2
#CREATED: 08 November 2021
#
#
#DESCRIPTION: This script install on your local server the latest stable
#   version of iSkyLIMS application 
#
#
#================================================================
# END_OF_HEADER
#================================================================

ISKYLIMS_VERSION="2.0.x"
. ./install_settings.txt

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

apache_check(){
    if ! pidof apache2 > /dev/null ; then
        # web server down, restart the server
        echo "Apache Server is down... Trying to restart Apache"
        systemctl restart apache2.service
        sleep 10
        if pidof apache2 > /dev/null ; then
            echo "Apache Server is up"
        else
            echo -e "${RED}ERROR : Unable to start Apache ${NC}"
            echo -e "${RED}ERROR : Solve the issue with Apache server and run again the installation script ${NC}"
            exit 1
        fi
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
apache_check
echo "Successful check for apache"
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
        libmysqlclient-dev apache2-dev vim libapache2-mod-wsgi-py3
    #apt-get install build-essential  -y 
    #apt-get install libghc-zlib-dev libbz2-dev libssl1.0-dev -y
    #apt-get install git libpango1.0 libpango1.0-dev   -y
fi

echo "Installing Interop"
mkdir -p /opt/interop
cd /opt/interop

wget https://github.com/Illumina/interop/releases/download/v1.1.15/InterOp-1.1.15-Linux-GNU.tar.gz
tar -xf  InterOp-1.1.15-Linux-GNU.tar.gz
ln -s InterOp-1.1.15-Linux-GNU interop
rm InterOp-1.1.15-Linux-GNU.tar.gz


echo "iSkyLIMS installation"
cd /opt/iSkyLIMS
git checkout master

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
cp conf/settings.py /opt/iSkyLIMS/iSkyLIMS/settings.py
cp conf/urls.py /opt/iSkyLIMS/iSkyLIMS/

sed -i "/^SECRET/c\\$(cat ~/.secret)" iSkyLIMS/settings.py
sed -i "s/djangouser/${DB_USER}/g" iSkyLIMS/settings.py
sed -i "s/djangopass/${DB_PASS}/g" iSkyLIMS/settings.py
sed -i "s/djangohost/${DB_SERVER_IP}/g" iSkyLIMS/settings.py
sed -i "s/djangoport/${DB_PORT}/g" iSkyLIMS/settings.py

sed -i "s/emailhost/${EMAIL_HOST}/g" iSkyLIMS/settings.py
sed -i "s/emailport/${EMAIL_PORT}/g" iSkyLIMS/settings.py
sed -i "s/emailhostuser/${EMAIL_HOST_USER}/g" iSkyLIMS/settings.py
sed -i "s/emailhostpassword/${EMAIL_HOST_PASSWORD}/g" iSkyLIMS/settings.py
sed -i "s/emailhosttls/${EMAIL_USE_TLS}/g" iSkyLIMS/settings.py
sed -i "s/localserverip/${LOCAL_SERVER_IP}/g" iSkyLIMS/settings.py


echo "Creating the database structure for iSkyLIMS"
python3 manage.py migrate
python3 manage.py makemigrations django_utils iSkyLIMS_core iSkyLIMS_wetlab iSkyLIMS_drylab iSkyLIMS_clinic
python3 manage.py migrate

python3 manage.py collectstatic

echo "Change owner of files to Apache user"
chown -R www-data:www-data /opt/iSkyLIMS

echo "Loading in database initial data"
python3 manage.py loaddata conf/new_installation_loading_tables.json

echo "Running crontab"
python3 manage.py crontab add
mv /var/spool/cron/crontabs/root /var/spool/cron/crontabs/www-data
chown www-data /var/spool/cron/crontabs/www-data 



echo "Updating Apache configuration"
cp conf/apache2.conf /etc/apache2/sites-available/000-default.conf
echo  'LoadModule wsgi_module "/opt/iSkyLIMS/virtualenv/lib/python3.8/site-packages/mod_wsgi/server/mod_wsgi-py38.cpython-38-x86_64-linux-gnu.so"' >/etc/apache2/mods-available/iskylims.load
cp conf/iskylims.conf /etc/apache2/mods-available/iskylims.conf

# Create needed symbolic links to enable the configurations:

ln -s /etc/apache2/mods-available/iskylims.load /etc/apache2/mods-enabled/
ln -s /etc/apache2/mods-available/iskylims.conf /etc/apache2/mods-enabled/

echo "Creating super user "
python3 manage.py createsuperuser

printf "\n\n%s"
printf "${BLUE}------------------${NC}\n"
printf "%s"
printf "${BLUE}Successfuly iSkyLIMS Installation version: ${ISKYLIMS_VERSION}${NC}\n"
printf "%s"
printf "${BLUE}------------------${NC}\n\n"




