#!/bin/bash

SCRIPT_VERSION=0.3

ISKYLIMS_VERSION="2.0.x"
. ./install_settings.txt

db_check(){
    mysqladmin -h $DB_SERVER_IP -u$DB_USER -p$DB_PASS -P$DB_PORT processlist > /dev/null ###user should have mysql permission on remote server.

    if ! [ $? -eq 0 ]; then
        echo -e "${RED}ERROR : Unable to connect to database. Check if your database is running and accessible${NC}"
        exit 1
    fi
    RESULT=`mysqlshow --user=$DB_USER --password=$DB_PASS --host=$DB_SERVER_IP --port=$DB_PORT | grep -o $DB_NAME`

    if  ! [ "$RESULT" == $DB_NAME ] ; then
        echo -e "${RED}ERROR : iSkyLIMS database is not defined yet ${NC}"
        echo -e "${RED}ERROR : Create iSkyLIMS database on your mysql server and run again the installation script ${NC}"
        exit 1
    fi
}

apache_check(){
    if [[ $linux_distribution == "Ubuntu" ]]; then
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
    elif [[ $linux_distribution == "CentOs" || $linux_distribution == "RedHatEnterprise" ]]; then
        if ! pidof httpd > /dev/null ; then
            # web server down, restart the server
            echo "Apache Server is down... Trying to restart Apache"
            systemctl restart httpd
            sleep 10
            if pidof httpd > /dev/null ; then
                echo "Apache Server is up"
            else
                echo -e "${RED}ERROR : Unable to start Apache ${NC}"
                echo -e "${RED}ERROR : Solve the issue with Apache server and run again the installation script ${NC}"
                exit 1
            fi
        fi
    fi
}
python_check(){

    python_version=$(su -c python3 --version $user)
    if [[ $python_version == "" ]]; then
        echo -e "${RED}ERROR : Python3 is not found in your system ${NC}"
        echo -e "${RED}ERROR : Solve the issue with Python and run again the installation script ${NC}"
        exit 1
    fi
    p_version=$(echo $python_version | cut -d"." -f2)
    if (( $p_version < 7 )); then
        echo -e "${RED}ERROR : Application requieres at least the version 3.7.x of Python3  ${NC}"
        echo -e "${RED}ERROR : Solve the issue with Python and run again the installation script ${NC}"
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

user=$SUDO_USER
group=$(groups | cut -d" " -f1)

#Linux distribution
linux_distribution=$(lsb_release -i | cut -f 2-)

echo "Checking main requirements"
python_check
printf "${BLUE}Valid version of Python${NC}\n"
db_check
printf "${BLUE}Successful check for database${NC}\n"
apache_check
printf "${BLUE}Successful check for apache${NC}\n"

#================================================================

read -p "Are you sure you want to install iSkyLIMS in this server? " -n 1 -r
echo    # (optional) move to a new line
if [[ ! $REPLY =~ ^[Yy]$ ]] ; then
    echo "Exiting without running iSkyLIMS installation"
    exit 1
fi

#================================================================
if [[ $linux_distribution == "Ubuntu" ]]; then
    echo "Software installation for Ubuntu"
    apt-get update && apt-get upgrade -y
    apt-get install -y \
        apt-utils wget \
        libmysqlclient-dev apache2-dev \
        python3-venv
fi

if [[ $linux_distribution == "CentOS" || $linux_distribution == "RedHatEnterprise" ]]; then
    echo "Software installation for Centos/RedHat"
    yum install zlib-devel bzip2-devel openssl-devel \
                wget httpd-devel mysql-libs
fi

echo "Installing Interop"
if [ -d /opt/interop ]; then
    echo "There is already an interop installation"
    echo "Skipping the Interop installation"
else
    mkdir -p /opt/interop
    cd /opt/interop
    echo "Downloading interop software"
    wget https://github.com/Illumina/interop/releases/download/v1.1.15/InterOp-1.1.15-Linux-GNU.tar.gz
    tar -xf  InterOp-1.1.15-Linux-GNU.tar.gz
    ln -s InterOp-1.1.15-Linux-GNU interop
    rm InterOp-1.1.15-Linux-GNU.tar.gz
    echo "Interop is now installed"
    cd -
fi

echo "Starting iSkyLIMS installation"
if [ -d $INSTALL_PATH ]; then
    echo "There already is an installation of relecov-platform in $INSTALL_PATH."
    read -p "Do you want to remove current installation and reinstall? (Y/N) " -n 1 -r
    echo    # (optional) move to a new line
    if [[ ! $REPLY =~ ^[Yy]$ ]] ; then
        echo "Exiting without running iSkyLIMS installation"
        exit 1
    else
        rm -rf $INSTALL_PATH
    fi
fi


## Create the installation folder
mkdir $INSTALL_PATH/

rsync -rlv README.md LICENSE conf iSkyLIMS_core iSkyLIMS_drylab \
        iSkyLIMS_wetlab iSkyLIMS_clinic django_utils $INSTALL_PATH

cd $INSTALL_PATH

## Create apache group if it does not exist.
if ! grep -q apache /etc/group
then
    groupadd apache
fi

## Fix permissions and owners

if [ $LOG_TYPE == "symbolic_link" ]; then
    if [ -d $LOG_PATH ]; then
        ln -s $LOG_PATH  $INSTALL_PATH/logs
    chmod 775 $LOG_PATH
    else
        echo "Log folder path: $LOG_PATH does not exist. Fix it in the install_settings.txt and run again."
    exit 1
    fi
else
    mkdir -p $INSTALL_PATH/logs
    chown $user:apache $INSTALL_PATH/logs
    chmod 775 $INSTALL_PATH/logs
fi

# Create necessary folders
mkdir -p $INSTALL_PATH/documents/wetlab/tmp
chown $user:apache $INSTALL_PATH/documents
chmod 775 $INSTALL_PATH/documents
chown $user:apache $INSTALL_PATH/documents/tmp
chmod 775 $INSTALL_PATH/documents
mkdir -p $INSTALL_PATH/documents/drylab
chown $user:apache $INSTALL_PATH/documents/drylab
chmod 775 $INSTALL_PATH/documents/drylab
echo "Created folders for logs and documents "


# install virtual environment
echo "Creating virtual environment"
if [ -d $INSTALL_PATH/virtualenv ]; then
    echo "There already is a virtualenv for iSkyLIMS in $INSTALL_PATH."
    read -p "Do you want to remove current virtualenv and reinstall? (Y/N) " -n 1 -r
    echo    # (optional) move to a new line
    if [[ ! $REPLY =~ ^[Yy]$ ]] ; then
        rm -rf $INSTALL_PATH/virtualenv
        bash -c "$PYTHON_BIN_PATH -m venv virtualenv"
    else
        echo "virtualenv alredy defined. Skipping."
    fi
else
    bash -c "$PYTHON_BIN_PATH -m venv virtualenv"
fi

echo "activate the virtualenv"
source virtualenv/bin/activate

# Starting iSkyLIMS
python3 -m pip install -r conf/pythonPackagesRequired.txt
django-admin startproject iSkyLIMS .
grep ^SECRET iSkyLIMS/settings.py > ~/.secret


# Copying config files and script
cp conf/settings.py $INSTALL_PATH/iSkyLIMS/settings.py
cp conf/urls.py $INSTALL_PATH/iSkyLIMS/

sed -i "/^SECRET/c\\$(cat ~/.secret)" $INSTALL_PATH/iSkyLIMS/settings.py
sed -i "s/djangouser/${DB_USER}/g" $INSTALL_PATH/iSkyLIMS/settings.py
sed -i "s/djangopass/${DB_PASS}/g" $INSTALL_PATH/iSkyLIMS/settings.py
sed -i "s/djangohost/${DB_SERVER_IP}/g" $INSTALL_PATH/iSkyLIMS/settings.py
sed -i "s/djangoport/${DB_PORT}/g" $INSTALL_PATH/iSkyLIMS/settings.py
sed -i "s/djangodbname/${DB_NAME}/g" $INSTALL_PATH/iSkyLIMS/settings.py

sed -i "s/emailhostserver/${EMAIL_HOST_SERVER}/g" $INSTALL_PATH/iSkyLIMS/settings.py
sed -i "s/emailport/${EMAIL_PORT}/g" $INSTALL_PATH/iSkyLIMS/settings.py
sed -i "s/emailhostuser/${EMAIL_HOST_USER}/g" $INSTALL_PATH/iSkyLIMS/settings.py
sed -i "s/emailhostpassword/${EMAIL_HOST_PASSWORD}/g" $INSTALL_PATH/iSkyLIMS/settings.py
sed -i "s/emailhosttls/${EMAIL_USE_TLS}/g" $INSTALL_PATH/iSkyLIMS/settings.py
sed -i "s/localserverip/${LOCAL_SERVER_IP}/g" $INSTALL_PATH/iSkyLIMS/settings.py


echo "Creating the database structure for iSkyLIMS"
python3 manage.py migrate
python3 manage.py makemigrations django_utils iSkyLIMS_core iSkyLIMS_wetlab iSkyLIMS_drylab iSkyLIMS_clinic
python3 manage.py migrate

python3 manage.py collectstatic

# echo "Change owner of files to Apache user"
# chown -R www-data:www-data /opt/iSkyLIMS

echo "Loading in database initial data"
python3 manage.py loaddata conf/new_installation_loading_tables.json

echo "Running crontab"
python3 manage.py crontab add
mv /var/spool/cron/crontabs/root /var/spool/cron/crontabs/www-data
chown www-data /var/spool/cron/crontabs/www-data



echo "Updating Apache configuration"
if [[ $linux_distribution == "Ubuntu" ]]; then
    cp conf/apache2.conf /etc/apache2/sites-available/000-default.conf
    echo  'LoadModule wsgi_module "/opt/iSkyLIMS/virtualenv/lib/python3.10.6/site-packages/mod_wsgi/server/mod_wsgi-py38.cpython-38-x86_64-linux-gnu.so"' >/etc/apache2/mods-available/iskylims.load
    cp conf/iskylims.conf /etc/apache2/mods-available/iskylims.conf

    # Create needed symbolic links to enable the configurations:

    ln -s /etc/apache2/mods-available/iskylims.load /etc/apache2/mods-enabled/
    ln -s /etc/apache2/mods-available/iskylims.conf /etc/apache2/mods-enabled/
fi

if [[ $linux_distribution == "CentOS" || $linux_distribution == "RedHatEnterprise" ]]; then
    cp conf/iskylims_apache_centos_redhat.conf /etc/httpd/conf.d/iskylims.conf
fi


echo "Creating super user "
echo "User name must be admin"
python3 manage.py createsuperuser

printf "\n\n%s"
printf "${BLUE}------------------${NC}\n"
printf "%s"
printf "${BLUE}Successfuly iSkyLIMS Installation version: ${ISKYLIMS_VERSION}${NC}\n"
printf "%s"
printf "${BLUE}------------------${NC}\n\n"




