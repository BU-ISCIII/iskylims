#!/bin/bash

ISKYLIMS_VERSION="2.x.x"

usage() {
cat << EOF
This script install and upgrade the iskylims app.

usage : $0 --upgrade --dev --conf
    Optional input data:
    --install | Install iskylims full/dep/app
    --upgrade | Upgrade iskylims full/dep/app
    --dev     | Use the develop version instead of main release
    --conf    | Select custom configuration file. Default: ./install_settings.txt
    --tables  | Load the first inital tables for upgrades in conf folder
    --script  | Run a migration script.
    --ren_app | Rename apps required for the upgrade migration to 2.3.1 


Examples:
    Install iskylims only dep
    sudo $0 --install dep

    Install only iSkyLIMS app
    $0 --install app

    Upgrade using develop code
    $0 --upgrade full --dev

    Upgrade running migration script and update initial tables
    $0 --upgrade full --script <migration_script> --tables

    Make adjustments for apps renaming in upgrade 2.3.0 to 2.3.1
    $0 --upgrade full --ren_app --script <migration_script> --tables
EOF
}

db_check(){
    # user should have mysql permission on remote server.
    mysqladmin -h $DB_SERVER_IP -u$DB_USER -p$DB_PASS -P$DB_PORT processlist > /dev/null

    if ! [ $? -eq 0 ]; then
        echo -e "${RED}ERROR : Unable to connect to database. Check if your database is running and accessible${NC}"
        exit 1
    fi
    RESULT=`mysqlshow --user=$DB_USER --password=$DB_PASS --host=$DB_SERVER_IP --port=$DB_PORT | grep -o $DB_NAME`

    if  ! [ "$RESULT" == "$DB_NAME" ] ; then
        echo -e "${RED}ERROR : iskylims database is not defined yet ${NC}"
        echo -e "${RED}ERROR : Create iskylims database on your mysql server and run again the installation script ${NC}"
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

    python_version=$(su -c $PYTHON_BIN_PATH --version $user)
    if [[ $python_version == "" ]]; then
        echo -e "${RED}ERROR : Python3 is not found in your system ${NC}"
        echo -e "${RED}ERROR : Solve the issue with Python and run again the installation script ${NC}"
        exit 1
    fi
    p_version=$(echo $python_version | cut -d"." -f2)
    if (( $p_version < 7 )); then
        echo -e "${RED}ERROR : Application requieres at least the version 3.7.x of Python3  ${NC}"
        echo -e "${RED}ERROR : Solve the issue with python and run again the installation script ${NC}"
        exit 1
    fi
}

root_check(){
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
}

update_settings_and_urls(){
    # save SECRET KEY at home user directory
    grep ^SECRET $INSTALL_PATH/iskylims/settings.py > ~/.secret

    # Copying config files and script. TODO CHANGE iSkyLIMS to app name
    cp conf/template_settings.py $INSTALL_PATH/iskylims/settings.py
    cp conf/urls.py $INSTALL_PATH/iskylims
    
    # replacing dummy variables with real values
    sed -i "/^SECRET/c\\$(cat ~/.secret)" $INSTALL_PATH/iskylims/settings.py
    sed -i "s/djangouser/${DB_USER}/g" $INSTALL_PATH/iskylims/settings.py
    sed -i "s/djangopass/${DB_PASS}/g" $INSTALL_PATH/iskylims/settings.py
    sed -i "s/djangohost/${DB_SERVER_IP}/g" $INSTALL_PATH/iskylims/settings.py
    sed -i "s/djangoport/${DB_PORT}/g" $INSTALL_PATH/iskylims/settings.py
    sed -i "s/djangodbname/${DB_NAME}/g" $INSTALL_PATH/iskylims/settings.py

    sed -i "s/emailhostserver/${EMAIL_HOST_SERVER}/g" $INSTALL_PATH/iskylims/settings.py
    sed -i "s/emailport/${EMAIL_PORT}/g" $INSTALL_PATH/iskylims/settings.py
    sed -i "s/emailhostuser/${EMAIL_HOST_USER}/g" $INSTALL_PATH/iskylims/settings.py
    sed -i "s/emailhostpassword/${EMAIL_HOST_PASSWORD}/g" $INSTALL_PATH/iskylims/settings.py
    sed -i "s/emailhosttls/${EMAIL_USE_TLS}/g" $INSTALL_PATH/iskylims/settings.py
    sed -i "s/localserverip/${LOCAL_SERVER_IP}/g" $INSTALL_PATH/iskylims/settings.py
    sed -i "s/localhost/${DNS_URL}/g" $INSTALL_PATH/iskylims/settings.py    
}

upgrade_venv(){
    rsync -rlv conf/requirements.txt $INSTALL_PATH/conf/requirements.txt
    echo "activate the virtualenv"
    source virtualenv/bin/activate
    echo "Installing required python packages"
    python -m pip install --upgrade pip
    python -m pip install -r conf/requirements.txt
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
        --install)  set -- "$@" -i ;;
        --upgrade)  set -- "$@" -u ;;
        --script)   set -- "$@" -s ;;
        --tables)   set -- "$@" -t ;;
        --dev)      set -- "$@" -d ;;
        --conf)     set -- "$@" -c ;;
        --ren_app)  set -- "$@" -r ;;

    # ADITIONAL
        --help)     set -- "$@" -h ;;
        --version)  set -- "$@" -v ;;
    # PASSING VALUE IN PARAMETER
        *)          set -- "$@" "$arg" ;;
    esac
done

# SETTING DEFAULT VALUES
ren_app=false
tables=false
git_branch="main"
conf="./install_settings.txt"
install=true
install_type="full"
upgrade=false
upgrade_type="full"

# PARSE VARIABLE ARGUMENTS WITH getops
options=":c:s:i:u:drtvh"
while getopts $options opt; do
    case $opt in
        i ) 
            install=true
            upgrade=false
            if [[ "$OPTARG" -eq "full" || "$OPTARG" -eq "dep" || "$OPTARG" -eq "app" ]]; then
                install_type=$OPTARG
                upgrade_type=$OPTARG
            else
                echo "Upgrade is not set to one valid option. Use: --upgrade full/app/dep"
                exit 1
            fi
            ;;
        u )
            install=false
            upgrade=true
            if [[ "$OPTARG" -eq "full" || "$OPTARG" -eq "dep" || "$OPTARG" -eq "app" ]]; then
                upgrade_type=$OPTARG
                install_type=$OPTARG
            else
                echo "Upgrade is not set to one valid option. Use: --upgrade full/app/dep"
                exit 1
            fi
            ;;
        s )
            run_script=true
            migration_script=$OPTARG
            ;;
        t )
            tables=true
            ;;
        r )
            ren_app=true
            ;;
        d )
            git_branch="develop"
            ;;
        c )
            conf=$OPTARG
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
#=============================================================================
#                     SETTINGS CHECKINGS
#=============================================================================

if [ ! -f "$conf" ]; then
    printf "\n\n%s"
    printf "${RED}------------------${NC}\n"
    printf "${RED}Unable to start.${NC}\n"
    printf "${RED}Configuration File $conf does not exist.${NC}\n"
    printf "${RED}------------------${NC}\n"
    exit 1
fi

# Read configuration file

. $conf

# check if branch master/develop is defined and checkout
if [ "`git branch --list $git_branch`" ]; then
    git checkout $git_branch
else
    printf "\n\n%s"
    printf "${RED}------------------${NC}\n"
    printf "${RED}Unable to start.${NC}\n"
    printf "${RED}Git branch $git_branch is not define in ${PWD}.${NC}\n"
    printf "${RED}------------------${NC}\n"
    exit 1
fi

#================================================================
# CHECK REQUIREMENTS BEFORE STARTING INSTALLATION
#================================================================

echo "Checking main requirements"
python_check
printf "${BLUE}Valid version of Python${NC}\n"
db_check
printf "${BLUE}Successful check for database${NC}\n"
apache_check
printf "${BLUE}Successful check for apache${NC}\n"

if [ "$install_type" == "full" ] || [ "$install_type" == "dep" ] || [ "$upgrade_type" == "full" ] || [ "$upgrade_type" == "dep" ]; then
    printf "${YELLOW} Checking requirement of root  user when installation is full or dep ${NC}\n"
    root_check
    printf "${BLUE}Successful checking of root user${NC}\n"
fi

#=============================================================================
#                   UPGRADE INSTALLATION
# Check if parameter is passing to script to upgrade the installation
# If "upgrade" parameter is set then the script only execute the upgrade part.
# If other parameter as upgrade is given return usage message and exit
#=============================================================================

if [ $upgrade == true ]; then
    # check if upgrade keyword is given
    if [ ! -d $INSTALL_PATH ]; then
        printf "\n\n%s"
        printf "${RED}------------------${NC}\n"
        printf "${RED}Unable to start the upgrade.${NC}\n"
        printf "${RED}Folder $INSTALL_PATH does not exist.${NC}\n"
        printf "${RED}------------------${NC}\n"
        exit 1
    fi
    #================================================================
    # MAIN_BODY FOR UPGRADE
    #================================================================
    printf "\n\n%s"
    printf "${YELLOW}------------------${NC}\n"
    printf "%s"
    printf "${YELLOW}Starting iSkyLIMS Upgrade version: ${ISKYLIMS_VERSION}${NC}\n"
    printf "%s"
    printf "${YELLOW}------------------${NC}\n\n"
    
    if [ "$upgrade_type" = "full" ] || [ "$upgrade_type" = "dep" ]; then
        if [ -d $INSTALL_PATH/virtualenv ]; then
            read -p "Do you want to remove current virtualenv and reinstall? (Y/N) " -n 1 -r
            echo    # (optional) move to a new line
            if [[ $REPLY =~ ^[Yy]$ ]] ; then
                rm -rf $INSTALL_PATH/virtualenv
                cd $INSTALL_PATH
                bash -c "$PYTHON_BIN_PATH -m venv virtualenv"
                upgrade_venv
                cd -
            else
                cd $INSTALL_PATH
                upgrade_venv
                cd -
            fi    
        else
            echo "There is no virtualenv to upgrade in $INSTALL_PATH."
            read -p "Do you want to create a new virtualenv and reinstall? (Y/N) " -n 1 -r
            echo    # (optional) move to a new line
            if [[ $REPLY =~ ^[Yy]$ ]] ; then
                cd $INSTALL_PATH
                bash -c "$PYTHON_BIN_PATH -m venv virtualenv"
                upgrade_venv
                cd -
            else
                echo "Exiting..."
                exit 0
            fi
        fi
    fi

    if [ "$upgrade_type" = "full" ] || [ "$upgrade_type" = "app" ]; then
        
        # Delete git and no copy files stuff
        if [ $ren_app == true ] ; then
            # remove all previous migrations and make a fake initial
            # delete existing migrations file
            rm -rf $INSTALL_PATH/django_utils/migrations/*
            rm -rf $INSTALL_PATH/iSkyLIMS_core/migrations/*
            rm -rf $INSTALL_PATH/iSkyLIMS_wetlab/migrations/*
            rm -rf $INSTALL_PATH/iSkyLIMS_drylab/migrations/*

            cd $INSTALL_PATH
            echo "activate the virtualenv"
            source virtualenv/bin/activate

            echo "Create a fake initial"
            python manage.py makemigrations django_utils iSkyLIMS_core \
                iSkyLIMS_wetlab iSkyLIMS_drylab
            python manage.py migrate --fake-initial

            if [ -d "$INSTALL_PATH/iSkyLIMS_core" ]; then
                echo "Changing app dir names in $INSTALL_PATH..."
                rm -rf $INSTALL_PATH/.git $INSTALL_PATH/.github $INSTALL_PATH/.gitignore \
                    $INSTALL_PATH/.Rhistory $INSTALL_PATH/docker-compose.yml $INSTALL_PATH/docker_iskylims_install.sh \
                    $INSTALL_PATH/Dockerfile $INSTALL_PATH/install.sh $INSTALL_PATH/install_settings.txt 
                mv $INSTALL_PATH/iSkyLIMS_core $INSTALL_PATH/core 
                mv $INSTALL_PATH/iSkyLIMS_wetlab $INSTALL_PATH/wetlab
                mv $INSTALL_PATH/iSkyLIMS_drylab $INSTALL_PATH/drylab
                mv $INSTALL_PATH/iSkyLIMS_clinic $INSTALL_PATH/clinic
                echo "Done changing app dir names in $INSTALL_PATH..."
            fi
            if [ -d "iSkyLIMS" ]; then
                mv iSkyLIMS/ iskylims/
                sed -i "s/iSkyLIMS/iskylims/g" $INSTALL_PATH/iskylims/wsgi.py
                sed -i "s/iSkyLIMS/iskylims/g" $INSTALL_PATH/manage.py
            fi
            cd -
        fi

        # update installation by sinchronize folders
        echo "Copying files to installation folder"
        rsync -rlv conf/ $INSTALL_PATH/conf/
        rsync -rlv --fuzzy --delay-updates --delete-delay \
              --exclude "logs" --exclude "documents" --exclude "migrations" --exclude "__pycache__" \
              README.md LICENSE conf core drylab clinic wetlab django_utils $INSTALL_PATH
        
        # update the settings.py and the main urls
        echo "Update settings and url file."
        update_settings_and_urls
        # update illumina template files.# Copy illumina sample sheet templates
        mkdir -p $INSTALL_PATH/documents/wetlab/templates/
        cp $INSTALL_PATH/conf/*_template.csv $INSTALL_PATH/documents/wetlab/templates/

        cd $INSTALL_PATH
        echo "activate the virtualenv"
        source virtualenv/bin/activate
        ### RENAME APP  in database and migration files ####
        if [ $ren_app == true ] ; then
            
            echo "Modifying names in migration files..."
            sed -i 's/iSkyLIMS_core/core/g' core/migrations/*.py
            # sed -i 's/iSkyLIMS_clinic/clinic/g' clinic/migrations/*.py
            # sed -i 's/iSkyLIMS_core/core/g' clinic/migrations/*.py
            sed -i 's/iSkyLIMS_drylab/drylab/g' drylab/migrations/*.py
            sed -i 's/iSkyLIMS_wetlab/wetlab/g' drylab/migrations/*.py
            sed -i 's/iSkyLIMS_core/core/g' drylab/migrations/*.py
            sed -i 's/iSkyLIMS_wetlab/wetlab/g' wetlab/migrations/*.py
            sed -i 's/iSkyLIMS_drylab/drylab/g' wetlab/migrations/*.py
            sed -i 's/iSkyLIMS_core/core/g' wetlab/migrations/*.py
            sed -i 's/iSkyLIMS_core/core/g' core/migrations/*.py
            echo "Done modifying names in migration files..."
            
            echo "Modifying database names and constraints..."
            mysql -u $DB_USER -p$DB_PASS -D $DB_NAME -h $DB_SERVER_IP \
                -e 'UPDATE django_content_type SET app_label = REPLACE(app_label , "iSkyLIMS_core", "core") WHERE app_label like ("iSkyLIMS_%");'
            # mysql -u $DB_USER -p$DB_PASS -D $DB_NAME -h $DB_SERVER_IP  \
            #    -e 'UPDATE django_content_type SET app_label = REPLACE(app_label , "iSkyLIMS_clinic", "clinic") WHERE app_label like ("iSkyLIMS_%");'
            mysql -u $DB_USER -p$DB_PASS -D $DB_NAME -h $DB_SERVER_IP \
                -e 'UPDATE django_content_type SET app_label = REPLACE(app_label , "iSkyLIMS_wetlab", "wetlab") WHERE app_label like ("iSkyLIMS_%");'
            mysql -u $DB_USER -p$DB_PASS -D $DB_NAME -h $DB_SERVER_IP \
                -e 'UPDATE django_content_type SET app_label = REPLACE(app_label , "iSkyLIMS_drylab", "drylab") WHERE app_label like ("iSkyLIMS_%");'
            
            mysql -u $DB_USER -p$DB_PASS -D $DB_NAME -h $DB_SERVER_IP \
                -e 'UPDATE django_migrations SET app = REPLACE(app , "iSkyLIMS_core", "core") WHERE app like ("iSkyLIMS_%");'
            # mysql -u $DB_USER -p$DB_PASS -D $DB_NAME -h $DB_SERVER_IP \
            #    -e 'UPDATE django_migrations SET app = REPLACE(app , "iSkyLIMS_clinic", "clinic") WHERE app like ("iSkyLIMS_%");'
            mysql -u $DB_USER -p$DB_PASS -D $DB_NAME -h $DB_SERVER_IP \
                -e 'UPDATE django_migrations SET app = REPLACE(app , "iSkyLIMS_wetlab", "wetlab") WHERE app like ("iSkyLIMS_%");'
            mysql -u $DB_USER -p$DB_PASS -D $DB_NAME -h $DB_SERVER_IP \
                -e 'UPDATE django_migrations SET app = REPLACE(app , "iSkyLIMS_drylab", "drylab") WHERE app like ("iSkyLIMS_%");'
            echo "Renaming tables"
            query_rename_table="SELECT CONCAT('RENAME TABLE ', TABLE_SCHEMA, '.', TABLE_NAME, \
                                ' TO ', TABLE_SCHEMA, '.', REPLACE(TABLE_NAME, 'iSkyLIMS_', ''), ';') \
                                AS query FROM information_schema.tables WHERE TABLE_SCHEMA = \"$DB_NAME\" AND TABLE_NAME LIKE 'iSkyLIMS_%';"
            mysql -u $DB_USER -p$DB_PASS -h $DB_SERVER_IP -e "$query_rename_table" \
                | xargs -I % echo "mysql -u$DB_USER -p'$DB_PASS' -D $DB_NAME -h $DB_SERVER_IP -e \"% \" " | bash
            echo "Renaming index"
            query_rename_unique_indexes="SELECT CONCAT('ALTER TABLE ', rcu.TABLE_SCHEMA, '.', rcu.TABLE_NAME, \
                                 ' RENAME INDEX ', rcu.CONSTRAINT_NAME, \
                                 ' TO ', REPLACE(rcu.CONSTRAINT_NAME, 'iSkyLIMS_', ''), ';') \
                                 AS query FROM information_schema.key_column_usage rcu \
                                 JOIN information_schema.table_constraints tc \
                                 ON tc.CONSTRAINT_NAME = rcu.CONSTRAINT_NAME WHERE rcu.TABLE_SCHEMA = \"$DB_NAME\" \
                                 AND rcu.CONSTRAINT_NAME LIKE 'iSkyLIMS_%' AND tc.CONSTRAINT_TYPE = 'UNIQUE' \
                                 GROUP BY rcu.TABLE_SCHEMA, rcu.TABLE_NAME, rcu.CONSTRAINT_NAME, tc.CONSTRAINT_TYPE, \
                                 rcu.REFERENCED_TABLE_SCHEMA, rcu.REFERENCED_TABLE_NAME;"
            mysql -u $DB_USER -p$DB_PASS -h $DB_SERVER_IP -e "$query_rename_unique_indexes"  \
                | xargs -I % echo "mysql -u$DB_USER -p'$DB_PASS' -D $DB_NAME -h $DB_SERVER_IP -e \"% \" " | bash
            echo "Renaming constraints"
            query_rename_constraints="SELECT CONCAT('ALTER TABLE ', rcu.TABLE_SCHEMA, '.', rcu.TABLE_NAME, \
                    ' DROP FOREIGN KEY ' , rcu.CONSTRAINT_NAME, ';', \
                    ' ALTER TABLE ', rcu.TABLE_SCHEMA, '.', rcu.TABLE_NAME, \
                    ' ADD CONSTRAINT ', REPLACE(rcu.CONSTRAINT_NAME, 'iSkyLIMS_', ''), ' ', \
                    tc.CONSTRAINT_TYPE, ' (', GROUP_CONCAT(rcu.COLUMN_NAME ORDER BY rcu.ORDINAL_POSITION SEPARATOR ', '), ')', \
                    IF(tc.CONSTRAINT_TYPE = 'FOREIGN KEY', \
                    CONCAT(' REFERENCES ', rcu.REFERENCED_TABLE_SCHEMA, '.', REPLACE(rcu.REFERENCED_TABLE_NAME, 'iSkyLIMS_', ''), ' (', \
                            GROUP_CONCAT(rcu.REFERENCED_COLUMN_NAME ORDER BY rcu.ORDINAL_POSITION SEPARATOR ', '), ') ON DELETE ', rc.DELETE_RULE), \
                    ''), ';') AS query \
                    FROM information_schema.key_column_usage rcu \
                    LEFT JOIN information_schema.table_constraints tc ON rcu.CONSTRAINT_NAME = tc.CONSTRAINT_NAME \
                    LEFT JOIN information_schema.referential_constraints rc ON rcu.CONSTRAINT_NAME = rc.CONSTRAINT_NAME \
                    WHERE rcu.TABLE_SCHEMA = '$DB_NAME' AND rcu.CONSTRAINT_NAME LIKE 'iSkyLIMS_%' \
                    GROUP BY rcu.TABLE_SCHEMA, rcu.TABLE_NAME, rcu.CONSTRAINT_NAME, tc.CONSTRAINT_TYPE, rcu.REFERENCED_TABLE_SCHEMA, rcu.REFERENCED_TABLE_NAME, rc.DELETE_RULE;"
            mysql -u $DB_USER -p$DB_PASS -h $DB_SERVER_IP -e "$query_rename_constraints" | xargs -I % echo "mysql -u$DB_USER -p'$DB_PASS' -D $DB_NAME -h $DB_SERVER_IP -e \"% \" " | bash

            echo "Done modifying database names and constraints..." 
                   
            # copy modified migration files
            echo "Copying custom migration files from conf."
            cp $INSTALL_PATH/conf/0002_core_migration_v2.3.1.py $INSTALL_PATH/core/migrations/0002_migration_v2_3_1.py
            cp $INSTALL_PATH/conf/0002_drylab_migration_v2.3.1.py $INSTALL_PATH/drylab/migrations/0002_migration_v2_3_1.py
            cp $INSTALL_PATH/conf/0002_wetlab_migration_v2.3.1.py $INSTALL_PATH/wetlab/migrations/0002_migration_v2_3_1.py
            # cp conf/0002_clinic_migration_v2.3.1.py clinic/migrations/0002_migration_v2_3_1.py
            cp $INSTALL_PATH/conf/0002_django_utils_migration_v2.3.1.py $INSTALL_PATH/django_utils/migrations/0002_migration_v2_3_1.py

            read -p "Do you want to proceed with the migrate command? (Y/N) " -n 1 -r
            echo    # (optional) move to a new line
            if [[ ! $REPLY =~ ^[Yy]$ ]] ; then
                echo "Exiting without running migrate command."
                exit 1
            fi

            echo "activate the virtualenv"
            source virtualenv/bin/activate
            echo "Running migrate..."
            python manage.py migrate
            echo "Done migrate command."

        else
            echo "checking for database changes"
            if python manage.py makemigrations | grep -q "No changes"; then
                echo "No migration is required"
            else
                read -p "Do you want to proceed with the migrate command? (Y/N) " -n 1 -r
                echo    # (optional) move to a new line
                if [[ ! $REPLY =~ ^[Yy]$ ]] ; then
                    echo "Exiting without running migrate command."
                    exit 1
                fi
                echo "Running migrate..."
                python manage.py migrate
                echo "Done migrate command."
            fi
        fi     
        
        echo "Running collect statics..."
        python manage.py collectstatic
        echo "Done collect statics"
        
        if [ $tables == true ] ; then
            echo "Loading pre-filled tables..."
            python manage.py loaddata conf/first_install_tables.json
            echo "Done loading pre-filled tables..."
        fi

        if [ $run_script ]; then
            echo "Running migration script: $migration_script"
            python manage.py runscript $migration_script
            echo "Done migration script: $migration_script"
        fi

        cd -

        # Linux distribution
        linux_distribution=$(lsb_release -i | cut -f 2-)

        echo ""
        echo "Restart apache server to update changes"
        if [[ $linux_distribution == "Ubuntu" ]]; then
            apache_user="apache2"
        else
            apache_user="httpd"
        fi
        
        # systemctl restart $apache_user

        if ! [ $? -eq 0 ]; then
            echo -e "${ORANGE}Apache server restart failed. trying with sudo{NC}"
            sudo systemctl restart $apache_user
        fi
    fi
    printf "\n\n%s"
    printf "${BLUE}------------------${NC}\n"
    printf "%s"
    printf "${BLUE}Successfuly upgrade of iSKyLIMS version: ${ISKYLIMS_VERSION}${NC}\n"
    printf "%s"
    printf "${BLUE}------------------${NC}\n\n"    
    # exit once upgrade is finished
    exit 0

fi

#================================================================
# INSTALL REPOSITORY REQUIRED SOFTWARE AND PYTHON VIRTUAL ENVIRONMENT
#================================================================

if [ $install == true ]; then

    if [ "$install_type" == "full" ] || [ "$install_type" == "dep" ]; then

        #================================================================
        # MAIN_BODY FOR INSTALL
        #================================================================
        printf "\n\n%s"
        printf "${YELLOW}------------------${NC}\n"
        printf "%s"
        printf "${YELLOW}Starting iSkyLIMS install version: ${ISKYLIMS_VERSION}${NC}\n"
        printf "%s"
        printf "${YELLOW}------------------${NC}\n\n"

        user=$SUDO_USER
        group=$(groups | cut -d" " -f1)

        # Find out server Linux distribution
        linux_distribution=$(lsb_release -i | cut -f 2-)

        read -p "Are you sure you want to install software required dependences for iSkyLIMS? (Y/N) " -n 1 -r
            echo    # (optional) move to a new line
            if [[ ! $REPLY =~ ^[Yy]$ ]] ; then
                echo "Exiting without installing required dependences software for iSkyLIMS installation"
                exit 1
            fi

        echo "Installing Interop"
        if [ -d /opt/interop ]; then
            echo "There is already an interop installation"
            echo "Skipping Interop installation"
        else
            cd /opt
            echo "Downloading interop software"
            wget https://github.com/Illumina/interop/releases/download/v1.1.15/InterOp-1.1.15-Linux-GNU.tar.gz
            tar -xf  InterOp-1.1.15-Linux-GNU.tar.gz
            ln -s InterOp-1.1.15-Linux-GNU interop
            rm InterOp-1.1.15-Linux-GNU.tar.gz
            echo "Interop is now installed"
            cd -
        fi

        if [[ $linux_distribution == "Ubuntu" ]]; then
            echo "Software installation for Ubuntu"
            apt-get update && apt-get upgrade -y
            apt-get install -y \
                apt-utils wget \
                libmysqlclient-dev \
                python3-venv  \
                libpq-dev \
                python3-dev python3-pip python3-wheel \
                apache2-dev
        fi

        if [[ $linux_distribution == "CentOS" || $linux_distribution == "RedHatEnterprise" ]]; then
            echo "Software installation for Centos/RedHat"
            yum groupinstall "Development tools"
            yum install zlib-devel bzip2-devel openssl-devel \
                        wget httpd-devel mysql-libs sqlite sqlite-devel \
                        mariadb-devel libffi-devel
        fi

        # install virtual environment
        echo "Creating virtual environment"
        if [ -d $INSTALL_PATH/virtualenv ]; then
            echo "There already is a virtualenv for iskylims in $INSTALL_PATH."
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

        # Install python packages required for iSkyLIMS
        echo "Installing required python packages"
        python -m pip install wheel
        python -m pip install -r conf/requirements.txt

        ## Create apache group if it does not exist.
        if ! grep -q apache /etc/group
        then
            groupadd apache
        fi

        if [ install == "full" || install == "app" ]; then
            echo "Software dep are successfuly installed"
        else
            printf "\n\n%s"
            printf "${BLUE}------------------${NC}\n"
            printf "%s"
            printf "${BLUE}Software dep are successfuly installed${NC}\n"
            printf "%s"
            printf "${BLUE}------------------${NC}\n\n"
            exit 0
        fi
    fi

    #================================================================
    # INSTALL iSkyLIMS PLATFORM APPLICATION
    #================================================================

    if [ "$install_type" == "full" ] || [ "$install_type" == "app" ]; then

        read -p "Are you sure you want to install iSkyLIMS app in this server? (Y/N) " -n 1 -r
        echo    # (optional) move to a new line
        if [[ ! $REPLY =~ ^[Yy]$ ]] ; then
            echo "Exiting without installing required software for iSkyLIMS installation"
            exit 1
        fi

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

        echo "Starting iSkyLIMS installation"
        if [ -d $INSTALL_PATH ]; then
            echo "There already is an installation of iskylims in $INSTALL_PATH."
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
        mkdir -p $INSTALL_PATH

        rsync -rlv README.md LICENSE conf core drylab \
                wetlab clinic django_utils $INSTALL_PATH

        cd $INSTALL_PATH

        # Create necessary folders
        echo "Created documents structure"
        mkdir -p $INSTALL_PATH/documents/wetlab/tmp
        mkdir -p $INSTALL_PATH/documents/wetlab/SampleSheets
        mkdir -p $INSTALL_PATH/documents/wetlab/images_plot
        mkdir -p $INSTALL_PATH/documents/wetlab/templates
        chown $user:apache $INSTALL_PATH/documents
        chmod 775 $INSTALL_PATH/documents
        chown $user:apache $INSTALL_PATH/documents/wetlab/tmp
        chmod 775 $INSTALL_PATH/documents/wetlab/tmp
        chown $user:apache $INSTALL_PATH/documents/wetlab/SampleSheets
        chmod 775 $INSTALL_PATH/documents/wetlab/SampleSheets
        chown $user:apache $INSTALL_PATH/documents/wetlab/images_plot
        chmod 775 $INSTALL_PATH/documents/wetlab/images_plot
        chown $user:apache $INSTALL_PATH/documents/wetlab/templates
        chmod 775 $INSTALL_PATH/documents/wetlab/templates
        mkdir -p $INSTALL_PATH/documents/drylab
        chown $user:apache $INSTALL_PATH/documents/drylab
        chmod 775 $INSTALL_PATH/documents/drylab

        # Copy illumina sample sheet templates
        cp $INSTALL_PATH/conf/*_template.csv $INSTALL_PATH/documents/wetlab/templates/

        # Starting iSkyLIMS
        echo "activate the virtualenv"
        source virtualenv/bin/activate

        echo "Creating iskylims project"
        django-admin startproject iskylims .
        
        # update the settings.py and the main urls
        update_settings_and_urls

        echo "Creating the database structure for iSkyLIMS"
        python manage.py migrate
        python manage.py makemigrations django_utils core wetlab drylab clinic
        python manage.py migrate

        # copy static files 
        echo "Run collectstatic"
        python manage.py collectstatic

        echo "Loading in database initial data"
        python manage.py loaddata conf/first_install_tables.json

        echo "Running crontab"
        ## TODO: CHECK THIS.
        python manage.py crontab add
        mv /var/spool/cron/crontabs/root /var/spool/cron/crontabs/www-data
        chown www-data /var/spool/cron/crontabs/www-data

        echo "Updating Apache configuration"
        if [[ $linux_distribution == "Ubuntu" ]]; then
            cp conf/iskylims_apache_ubuntu.conf /etc/apache2/sites-available/000-default.conf
        fi

        if [[ $linux_distribution == "CentOS" || $linux_distribution == "RedHatEnterprise" ]]; then
            cp conf/iskylims_apache_centos_redhat.conf /etc/httpd/conf.d/iskylims.conf
        fi

        echo "Creating super user "
        python manage.py createsuperuser --username admin

        printf "\n\n%s"
        printf "${BLUE}------------------${NC}\n"
        printf "%s"
        printf "${BLUE}Successfuly iSkyLIMS Installation version: ${ISKYLIMS_VERSION}${NC}\n"
        printf "%s"
        printf "${BLUE}------------------${NC}\n\n"

        echo "Installation completed"
        exit 0
    fi
fi

printf "\n\n%s"
printf "${RED}------------------${NC}\n"
printf "%s"
printf "${RED}Invalid installation parameters${NC}\n"
printf "%s"
printf "${RED}------------------${NC}\n\n"
echo "See the usage examples"
usage
exit 1
