#!/bin/bash

APP_VERSION="3.1.0dev"

# usage: prints the command line help and usage examples.
usage() {
cat << EOF
This script install and upgrade the iskylims app.

usage : $0 --upgrade --git_revision --conf
    Optional input data:
    --install       | Install iskylims full/dep/app
    --upgrade       | Upgrade iskylims full/dep/app
    --git_revision  | Git revision name to run (it can be git branch, git version tag or commit SHA)
    --conf          | Select custom configuration file. Default: ./install_settings.txt
    --tables        | Load the first inital tables (from conf folder)
    --skip_tables   | Skip loading initial tables (even during install)
    --script        | Run a migration script after migrations.
    --script_before | Run a migration script before migrations.
    --script_after  | Run a migration script after migrations (same as --script).
    --ren_app       | Rename apps required for the upgrade migration to 3.0.0
    --docker        | Deprecated. Use --skip_apache_restart to avoid Apache checks/restart.


Examples:
    Install iskylims only dep
    sudo $0 --install dep

    Install only iSkyLIMS app
    $0 --install app

    Upgrade using develop code
    $0 --upgrade full --git_revision develop

    Upgrade running migration script and update initial tables
    $0 --upgrade full --script <migration_script> --tables

    Make adjustments for apps renaming in upgrade 2.3.0 to 2.3.1
    $0 --upgrade full --ren_app --script <migration_script> --tables

    Upgrade running pre/post migration scripts:
    $0 --upgrade app --script_before <pre_script> --script_after <post_script>
EOF
}

# log: write timestamped log entries to stdout.
_log_compose_entry() {
    local level="$1"; shift
    local message="$*"
    local timestamp
    timestamp="$(date '+%Y-%m-%d %H:%M:%S')"
    printf "%s [%s] %s" "$timestamp" "$level" "$message"
}

log() {
    local level="$1"; shift
    local message="$*"
    local entry
    entry="$(_log_compose_entry "$level" "$message")"
    printf "%s\n" "$entry"
}

# db_check: verifies connectivity to the configured MySQL instance using mysqladmin/mysqlshow.
db_check(){
    log "INFO" "Checking database connectivity against $DB_SERVER_IP:$DB_PORT"
    local mysqladmin_bin
    local mysqlshow_bin
    mysqladmin_bin="$(command -v mysqladmin || command -v mariadb-admin || true)"
    mysqlshow_bin="$(command -v mysqlshow || command -v mariadb-show || true)"

    if [ -z "$mysqladmin_bin" ] || [ -z "$mysqlshow_bin" ]; then
        log "ERROR" "mysql client tools not found (mysqladmin/mysqlshow or mariadb-admin/mariadb-show)."
        exit 1
    fi

    "$mysqladmin_bin" -h $DB_SERVER_IP -u$DB_USER -p$DB_PASS -P$DB_PORT processlist > /dev/null

    if ! [ $? -eq 0 ]; then
        log "ERROR" "Unable to connect to database. Check if your database is running and accessible"
        exit 1
    fi
    RESULT=`"$mysqlshow_bin" --user=$DB_USER --password=$DB_PASS --host=$DB_SERVER_IP --port=$DB_PORT | grep -o $DB_NAME`

    if  ! [ "$RESULT" == "$DB_NAME" ] ; then
        log "ERROR" "iskylims database is not defined yet"
        log "ERROR" "Create iskylims database on your mysql server and run again the installation script"
        exit 1
    fi
}

# apache_check: ensures apache/httpd service is running depending on distribution.
apache_check(){
    if [[ $linux_distribution == "Ubuntu" ]]; then
        if ! pidof apache2 > /dev/null ; then
            log "WARN" "Apache Server is down... Trying to restart Apache"
            systemctl restart apache2.service
            sleep 10
            if pidof apache2 > /dev/null ; then
                log "INFO" "Apache Server is up"
            else
                log "ERROR" "Unable to start Apache"
                log "ERROR" "Solve the issue with Apache server and run again the installation script"
                exit 1
            fi
        fi
    elif [[ $linux_distribution == "CentOs" || $linux_distribution == "RedHatEnterprise" ]]; then
        if ! pidof httpd > /dev/null ; then
            log "WARN" "Apache Server is down... Trying to restart Apache"
            systemctl restart httpd
            sleep 10
            if pidof httpd > /dev/null ; then
                log "INFO" "Apache Server is up"
            else
                log "ERROR" "Unable to start Apache"
                log "ERROR" "Solve the issue with Apache server and run again the installation script"
                exit 1
            fi
        fi
    fi
}

# python_check: confirm required Python version is available in PYTHON_BIN_PATH.
python_check(){
    python_version=$(su -c $PYTHON_BIN_PATH --version $user 2>&1)
    if [[ $python_version == "" ]]; then
        log "ERROR" "Python3 is not found in your system"
        log "ERROR" "Solve the issue with Python and run again the installation script"
        exit 1
    fi
    p_version=$(echo $python_version | cut -d"." -f2)
    if (( $p_version < 7 )); then
        log "ERROR" "Application requires at least version 3.7.x of Python3"
        log "ERROR" "Solve the issue with python and run again the installation script"
        exit 1
    fi
}

# root_check: enforce running privileged sections as root.
root_check(){
    if [[ $EUID -ne 0 ]]; then
        log "ERROR" "Exiting installation. This script must be run as root"
        exit 1
    fi
}

# update_settings_and_urls: rewrite Django settings and urls with deployment values.
update_settings_and_urls(){
    log "INFO" "Updating settings.py and urls.py with deployment values"
    grep ^SECRET $INSTALL_PATH/iskylims/settings.py > ~/.secret

    cp conf/template_settings.txt $INSTALL_PATH/iskylims/settings.py
    cp conf/urls.py $INSTALL_PATH/iskylims
    
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

# restore_git_ref: reset repository to branch/tag/commit active before script ran.
restore_git_ref() {
    echo "Restoring to initial git reference: $initial_git_ref"
    git checkout "$initial_git_ref" --quiet
}

# load_tables: wrapper to call Django loaddata with optional verbosity.
load_tables() {
    # Function parameters
    local data_file="${1:-conf/first_install_tables.json}"
    local verbose="${2:-false}"

    # Check if the file exists
    if [[ ! -f "$data_file" ]]; then
        echo "Error: The data file '$data_file' does not exist."
        return 1
    fi

    # Conditional message based on verbose mode
    if [[ "$verbose" == true ]]; then
        echo "Loading pre-filled tables from file: $data_file"
    fi

    # Load pre-filled tables
    python manage.py loaddata "$data_file"
    if [[ $? -eq 0 ]]; then
        echo "Tables loaded successfully from '$data_file'."
    else
        echo "Error loading tables from '$data_file'."
        return 1
    fi

    if [[ "$verbose" == true ]]; then
        echo "Table loading process completed."
    fi
}

# Ensure to recover current git branch/tag/SHA on script exit
initial_git_ref=$(git rev-parse --abbrev-ref HEAD || git rev-parse HEAD)
trap restore_git_ref EXIT

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
ORANGE='\033[0;33m'

# log_section: print a visually separated header in both console and log.
log_section() {
    local message="$1"
    log "INFO" "$message"
    printf "\n\n%s\n" "${YELLOW}------------------${NC}"
    printf "%b\n" "${YELLOW}${message}${NC}"
    printf "%s\n\n" "${YELLOW}------------------${NC}"
}

# log_info: convenience helper for blue info messages (console only).
log_info() {
    printf "%b\n" "${BLUE}$(_log_compose_entry "INFO" "$1")${NC}"
}

# log_warn: emit warning text in cyan for terminal visibility.
log_warn() {
    printf "%b\n" "${CYAN}$(_log_compose_entry "WARN" "$1")${NC}"
}

# log_error: emit error text in red for terminal visibility.
log_error() {
    printf "%b\n" "${RED}$(_log_compose_entry "ERROR" "$1")${NC}"
}

# abort_install: log an error and exit with optional status.
abort_install() {
    log_error "$1"
    exit "${2:-1}"
}

ensure_file_exists() {
    local file_path="$1"
    local friendly_name="${2:-$1}"
    if [ ! -f "$file_path" ]; then
        abort_install "Required file '$friendly_name' not found."
    fi
}

# load_install_config: source the selected install_settings file.
load_install_config() {
    ensure_file_exists "$conf" "$conf"
    # shellcheck disable=SC1090
    . "$conf"
}

# checkout_git_revision: ensure desired git revision exists and check it out safely.
checkout_git_revision() {
    if git rev-parse --verify "$git_branch" >/dev/null 2>&1; then
        if [[ $git_branch != $initial_git_ref ]]; then
            local local_changes
            local_changes=$(git status --porcelain)
            if [[ -n $local_changes ]]; then
                abort_install "Unable to switch to $git_branch. Commit or stash local changes first."
            fi
            printf "${YELLOW}Switching to revision %s.${NC}\n" "$git_branch"
            git checkout "$git_branch" --quiet
        else
            printf "${YELLOW}Using current revision: '%s'.${NC}\n" "$git_branch"
        fi
    else
        abort_install "Git reference $git_branch is not defined in ${PWD}."
    fi
}

# check_requirements: run Python/DB/Apache/root validations before install/upgrade.
check_requirements() {
    log_section "Checking main requirements"
    python_check
    log_info "Valid version of Python"
    if [[ "$operation_scope" == "full" || "$operation_scope" == "app" ]]; then
        db_check
        log_info "Successful check for database"
        if [ "$restart_apache" = true ]; then
            apache_check
            log_info "Successful check for apache"
        fi
    fi

    if [ "$install_type" == "full" ] || [ "$install_type" == "dep" ] || [ "$upgrade_type" == "full" ] || [ "$upgrade_type" == "dep" ]; then
        log_warn "Checking requirement of root user when installation is full or dep"
        root_check
        log_info "Successful checking of root user"
    fi
}

# rename_apps_if_needed: handles legacy app renaming and DB/migration adjustments when --ren_app is provided.
rename_apps_if_needed() {
    if [ $ren_app != true ]; then
        return 0
    fi

    rm -rf $INSTALL_PATH/django_utils/migrations/*
    rm -rf $INSTALL_PATH/iSkyLIMS_core/migrations/*
    rm -rf $INSTALL_PATH/iSkyLIMS_wetlab/migrations/*
    rm -rf $INSTALL_PATH/iSkyLIMS_drylab/migrations/*

    cd $INSTALL_PATH
    sed -i "s/ugettext/gettext/g" iSkyLIMS_wetlab/models.py
    sed -i "s/ugettext/gettext/g" iSkyLIMS_core/forms.py
    sed -i "s/ugettext/gettext/g" django_utils/forms.py
    echo "activate the virtualenv"
    source virtualenv/bin/activate

    echo "Create a fake initial"
    python manage.py makemigrations $FAKEINITIAL_MODULES
    python manage.py migrate --fake-initial

    if [ -d "$INSTALL_PATH/iSkyLIMS_core" ]; then
        echo "Changing app dir names in $INSTALL_PATH..."
        rm -rf $INSTALL_PATH/.git $INSTALL_PATH/.github $INSTALL_PATH/.gitignore \
            $INSTALL_PATH/.Rhistory $INSTALL_PATH/docker-compose.test.yml $INSTALL_PATH/docker_iskylims_install.sh \
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

    echo "Modifying database names and constraints..."
    mysql -u $DB_USER -p$DB_PASS -D $DB_NAME -h $DB_SERVER_IP \
        -e 'UPDATE django_content_type SET app_label = REPLACE(app_label , "iSkyLIMS_core", "core") WHERE app_label like ("iSkyLIMS_%");'
    mysql -u $DB_USER -p$DB_PASS -D $DB_NAME -h $DB_SERVER_IP \
        -e 'UPDATE django_content_type SET app_label = REPLACE(app_label , "iSkyLIMS_wetlab", "wetlab") WHERE app_label like ("iSkyLIMS_%");'
    mysql -u $DB_USER -p$DB_PASS -D $DB_NAME -h $DB_SERVER_IP \
        -e 'UPDATE django_content_type SET app_label = REPLACE(app_label , "iSkyLIMS_drylab", "drylab") WHERE app_label like ("iSkyLIMS_%");'

    mysql -u $DB_USER -p$DB_PASS -D $DB_NAME -h $DB_SERVER_IP \
        -e 'UPDATE django_migrations SET app = REPLACE(app , "iSkyLIMS_core", "core") WHERE app like ("iSkyLIMS_%");'
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

    echo "Modifying names in migration files..."
    sed -i 's/iSkyLIMS_core/core/g' */migrations/*.py
    sed -i 's/iSkyLIMS_drylab/drylab/g' */migrations/*.py
    sed -i 's/iSkyLIMS_wetlab/wetlab/g' */migrations/*.py
    echo "Done modifying names in migration files..."

    echo "Copying custom migration files from conf."
    cp $INSTALL_PATH/conf/0002_core_migration_v3.0.0.py $INSTALL_PATH/core/migrations/0002_migration_v3_0_0.py
    cp $INSTALL_PATH/conf/0002_drylab_migration_v3.0.0.py $INSTALL_PATH/drylab/migrations/0002_migration_v3_0_0.py
    cp $INSTALL_PATH/conf/0002_wetlab_migration_v3.0.0.py $INSTALL_PATH/wetlab/migrations/0002_migration_v3_0_0.py
    cp $INSTALL_PATH/conf/0002_django_utils_migration_v3.0.0.py $INSTALL_PATH/django_utils/migrations/0002_migration_v3_0_0.py

    read -p "Do you want to proceed with the migrate command? (Y/N) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]] ; then
        log "WARN" "Exiting without running migrate command."
        exit 1
    fi

    echo "activate the virtualenv"
    source virtualenv/bin/activate
    echo "Running migrate..."
    python manage.py migrate
    echo "Done migrate command."

    cd -
}

# install_system_packages: install InterOp and distro-specific OS packages required by iSkyLIMS.
install_system_packages() {
    if [ "${SKIP_SYSTEM_PACKAGES:-}" = "1" ]; then
        echo "Skipping system package installation (SKIP_SYSTEM_PACKAGES=1)"
        return
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

    linux_distribution=$(lsb_release -i | cut -f 2-)

    if [[ $linux_distribution == "Ubuntu" ]]; then
        echo "Software installation for Ubuntu"
        apt-get update && apt-get upgrade -y
        apt-get install -y \
            apt-utils wget \
            libmysqlclient-dev \
            python3-venv  \
            libpq-dev \
            python3-dev python3-pip python3-wheel \
            apache2-dev cifs-utils \
            gnuplot

    elif [[ $linux_distribution == "CentOS" || $linux_distribution == "RedHatEnterprise" ]]; then
        echo "Software installation for Centos/RedHat"
        yum groupinstall "Development tools"
        yum install zlib-devel bzip2-devel openssl-devel \
                    wget httpd-devel mysql-libs sqlite sqlite-devel \
                    mariadb-devel libffi-devel \
                    gnuplot cifs-utils
    fi
}

# run_django_deploy: execute makemigrations/migrate and optional fixture/superuser steps.
run_django_deploy() {
    local mode="${1:-install}"
    if [ "$run_script_before" = true ]; then
        for val in "${migration_script_before[@]}"; do
            if [[ $val = *","* ]]; then
                parameters=(${val//,/ })
                echo "Running pre-migration script: ${parameters[0]}"
                ./manage.py runscript ${parameters[0]} --script-args ${parameters[1]}
                echo "Done pre-migration script: ${parameters[0]}"
            else
                echo "Running pre-migration script: $val"
                ./manage.py runscript $val
                echo "Done pre-migration script: $val"
            fi
        done
    fi

    if [ "$mode" = "upgrade" ]; then
        echo "Applying migrations in fake-initial mode"
        python manage.py migrate --noinput --fake-initial
        # Second pass ensures non-initial migrations are applied after fake-initial.
        echo "Applying migrations"
        python manage.py migrate --noinput
    else
        echo "Applying migrations"
        python manage.py migrate --noinput
    fi

    if [ "$tables" = true ]; then
        echo "Loading pre-filled tables..."
        load_tables "$prefilled_tables" true
        echo "Done loading pre-filled tables..."
    fi

    if [ "$run_script" = true ]; then
        for val in "${migration_script[@]}"; do
            if [[ $val = *","* ]]; then
                parameters=(${val//,/ })
                echo "Running post-migration script: ${parameters[0]}"
                ./manage.py runscript ${parameters[0]} --script-args ${parameters[1]}
                echo "Done post-migration script: ${parameters[0]}"
            else
                echo "Running post-migration script: $val"
                ./manage.py runscript $val
                echo "Done post-migration script: $val"
            fi
        done
    fi

    if [ "$mode" = "install" ]; then
        echo "Creating super user "
        python manage.py createsuperuser --username admin
    fi
}

# sync_requirements_file: copy repository requirements into the target installation path.
sync_requirements_file() {
    mkdir -p $INSTALL_PATH/conf
    rsync -rlv conf/requirements.txt $INSTALL_PATH/conf/requirements.txt
}

# setup_virtualenv: create or refresh the Python virtualenv depending on mode.
setup_virtualenv() {
    local mode="$1"
    cd $INSTALL_PATH
    if [ "$mode" = "install" ]; then
        if [ -d virtualenv ]; then
            echo "There already is a virtualenv for iskylims in $INSTALL_PATH."
            read -p "Do you want to remove current virtualenv and reinstall? (Y/N) " -n 1 -r
            echo
            if [[ ! $REPLY =~ ^[Yy]$ ]] ; then
                rm -rf $INSTALL_PATH/virtualenv
                bash -c "$PYTHON_BIN_PATH -m venv virtualenv"
            else
                echo "virtualenv already defined. Skipping."
            fi
        else
            bash -c "$PYTHON_BIN_PATH -m venv virtualenv"
        fi
    else
        if [ -d virtualenv ]; then
            read -p "Do you want to remove current virtualenv and reinstall? (Y/N) " -n 1 -r
            echo
            if [[ $REPLY =~ ^[Yy]$ ]] ; then
                rm -rf $INSTALL_PATH/virtualenv
                bash -c "$PYTHON_BIN_PATH -m venv virtualenv"
            fi
        else
            read -p "There is no virtualenv. Do you want to create a new one? (Y/N) " -n 1 -r
            echo
            if [[ $REPLY =~ ^[Yy]$ ]] ; then
                bash -c "$PYTHON_BIN_PATH -m venv virtualenv"
            else
                echo "Exiting..."
                exit 0
            fi
        fi
    fi
    cd -
}

# prepare_documents_structure: ensure document directories and templates exist with correct permissions.
prepare_documents_structure() {
    echo "Created documents structure"
    mkdir -p $INSTALL_PATH/documents/wetlab
    mkdir -p $INSTALL_PATH/documents/wetlab/tmp
    mkdir -p $INSTALL_PATH/documents/wetlab/sample_sheet
    mkdir -p $INSTALL_PATH/documents/wetlab/images_plot
    mkdir -p $INSTALL_PATH/documents/wetlab/templates
    mkdir -p $INSTALL_PATH/documents/wetlab/sample_sheets_lib_prep
    mkdir -p $INSTALL_PATH/documents/drylab
    mkdir -p $INSTALL_PATH/documents/drylab/service_files

    chown -R $user:$apache_group $INSTALL_PATH/documents
    chmod 775 $INSTALL_PATH/documents

    cp $INSTALL_PATH/conf/*_template.csv $INSTALL_PATH/documents/wetlab/templates/
    cp $INSTALL_PATH/conf/samples_template.xlsx $INSTALL_PATH/documents/wetlab/templates/

    mkdir -p $INSTALL_PATH/documents/wetlab/collection_index_kits/
    cp $INSTALL_PATH/conf/collection_index_kits/*.txt $INSTALL_PATH/documents/wetlab/collection_index_kits/

    cp $INSTALL_PATH/conf/template_logging_config.ini $INSTALL_PATH/wetlab/logging_config.ini
    sed -i "s|INSTALL_PATH|${INSTALL_PATH}|g" $INSTALL_PATH/wetlab/logging_config.ini
}

# install_python_requirements: activate the venv and install required Python packages.
install_python_requirements() {
    cd $INSTALL_PATH
    echo "activate the virtualenv"
    source virtualenv/bin/activate
    echo "Installing required python packages"
    python -m pip install --upgrade pip
    python -m pip install wheel
    python -m pip install -r conf/requirements.txt
    cd -
}

ensure_virtualenv_ready() {
    if [ ! -d "$INSTALL_PATH/virtualenv" ]; then
        log_warn "Virtualenv missing. INSTALL_PATH=$INSTALL_PATH"
        ls -la "$INSTALL_PATH" || true
        abort_install "Virtualenv not found at $INSTALL_PATH/virtualenv. Run --install dep first."
    fi
}

# restart_apache_service: restart Apache/HTTPD unless running inside Docker or explicitly skipped.
restart_apache_service() {
    linux_distribution=$(lsb_release -i | cut -f 2-)
    if [[ $linux_distribution == "Ubuntu" ]]; then
        apache_daemon="apache2"
    else
        apache_daemon="httpd"
    fi
    if ! systemctl restart $apache_daemon; then
        echo -e "${ORANGE}Apache server restart failed. trying with sudo${NC}"
        sudo systemctl restart $apache_daemon
    fi
}

# run_dependency_stage: execute the dependency portion (system packages + venv + pip) for install or upgrade.
run_dependency_stage() {
    local mode="$1"

    if [ "$mode" = "install" ]; then
    log_section "Preparing dependency environment for installation"
        if [ -d $INSTALL_PATH ]; then
            echo "There already is an installation of iskylims in $INSTALL_PATH."
            read -p "Do you want to remove current installation and reinstall? (Y/N) " -n 1 -r
            echo
            if [[ ! $REPLY =~ ^[Yy]$ ]] ; then
                echo "Exiting without running iSkyLIMS installation"
                exit 1
            else
                rm -rf $INSTALL_PATH
            fi
        fi
        install_system_packages
        mkdir -p $INSTALL_PATH
        linux_distribution=$(lsb_release -i | cut -f 2-)
        if [[ $linux_distribution == "Ubuntu" ]]; then
            apache_group="www-data"
        else
            apache_group="apache"
        fi
        chown -R $user:$apache_group $INSTALL_PATH
        chmod 775 $INSTALL_PATH
    else
        log_section "Preparing dependency environment for upgrade"
        if [ ! -d $INSTALL_PATH ]; then
            abort_install "Unable to start the upgrade. Folder $INSTALL_PATH does not exist."
        fi
        install_system_packages
    fi

    sync_requirements_file
    setup_virtualenv "$mode"
    install_python_requirements
}

# upgrade_application_files: sync code/config and run upgrade-specific tasks (renames, migrations).
upgrade_application_files() {
    if [ ! -d $INSTALL_PATH ]; then
        abort_install "Unable to start the upgrade. Folder $INSTALL_PATH does not exist."
    fi

    log_section "Starting iSkyLIMS Upgrade version: ${APP_VERSION}"

    rename_apps_if_needed

    echo "Copying files to installation folder"
    rsync -rlv conf/ $INSTALL_PATH/conf/
    rsync -rlv --fuzzy --delay-updates --delete-delay \
          --exclude "logs" --exclude "documents" --exclude "__pycache__" \
          README.md LICENSE test conf $REQUIRED_MODULES $INSTALL_PATH

    cd $INSTALL_PATH
    ensure_virtualenv_ready
    echo "activate the virtualenv"
    source virtualenv/bin/activate

    if [ ! -f "$INSTALL_PATH/manage.py" ]; then
        echo "manage.py not found. Creating iskylims project"
        "$INSTALL_PATH/virtualenv/bin/python" -m django startproject iskylims .
    fi

    echo "Update settings and url file."
    update_settings_and_urls
    prepare_documents_structure

    run_django_deploy "upgrade"
    echo "Deleting static files..."
    if [ -d "$INSTALL_PATH/static" ]; then
        if command -v mountpoint >/dev/null 2>&1 && mountpoint -q "$INSTALL_PATH/static"; then
            echo "Static directory is a mount point. Skipping delete."
        else
            rm -rf "$INSTALL_PATH/static" || echo "Skipping static removal (busy)."
        fi
    fi
    echo "Running collect statics..."
    python manage.py collectstatic
    echo "Done collect statics"

    cd -
    log_section "Successfuly upgrade of iSKyLIMS version: ${APP_VERSION}"
}

# install_application_files: deploy Django project files, update settings, and run initial migrations.
install_application_files() {
    log_section "Starting iSkyLIMS install version: ${APP_VERSION}"

    user=${SUDO_USER:-$USER}
    group=$(groups | cut -d" " -f1)

    linux_distribution=$(lsb_release -i | cut -f 2-)

    if [[ $linux_distribution == "Ubuntu" ]]; then
        apache_group="www-data"
    else
        apache_group="apache"
    fi

    if [ "$install_type" == "full" ] || [ "$install_type" == "app" ]; then

        if [ $LOG_TYPE == "symbolic_link" ]; then
            if [ -d $LOG_PATH ]; then
                if [ -e "$INSTALL_PATH/logs" ]; then
                    echo "Log target $INSTALL_PATH/logs already exists. Leaving it unchanged."
                else
                    echo "Creating symbolic link to log folder"
                    ln -s "$LOG_PATH" "$INSTALL_PATH/logs"
                    chmod 775 "$LOG_PATH"
                fi
            else
                echo "Log folder path: $LOG_PATH does not exist. Fix it in the install_settings.txt and run again."
            exit 1
            fi
        else
            if  [ ! -d $INSTALL_PATH/logs ]; then
                mkdir -p $INSTALL_PATH/logs
                chown $user:$apache_group $INSTALL_PATH/logs
                chmod 775 $INSTALL_PATH/logs
            else
                echo "Log folder path: $INSTALL_PATH/logs already exist."
            fi
        fi

        rsync -rlv README.md LICENSE test conf $REQUIRED_MODULES $INSTALL_PATH

        cd $INSTALL_PATH

        prepare_documents_structure

        ensure_virtualenv_ready
        echo "activate the virtualenv"
        source virtualenv/bin/activate

        echo "Creating iskylims project"
        "$INSTALL_PATH/virtualenv/bin/python" -m django startproject iskylims .

        update_settings_and_urls

        run_django_deploy "install"

        echo "Run collectstatic"
        python manage.py collectstatic

        cd -

        log_section "Successfuly iSkyLIMS Installation version: ${APP_VERSION}"
        echo "Installation completed"
    fi
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
        --install)      set -- "$@" -i ;;
        --upgrade)      set -- "$@" -u ;;
        --script)       set -- "$@" -s ;;
        --script_before) set -- "$@" -p ;;
        --script_after) set -- "$@" -o ;;
        --script_prev)  set -- "$@" -p ;;
        --tables)       set -- "$@" -t ;;
        --skip_tables)  set -- "$@" -b ;;
        --git_revision) set -- "$@" -g ;;
        --conf)         set -- "$@" -c ;;
        --ren_app)      set -- "$@" -r ;;
        --docker)       set -- "$@" -k ;;
        --skip_apache_restart) set -- "$@" -a ;;

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
git_branch=$initial_git_ref
conf="./install_settings.txt"
install=true
install_type="full"
upgrade=false
upgrade_type="full"
docker=false
prefilled_tables="conf/first_install_tables.json"
restart_apache=true
run_script=false
run_script_before=false
migration_script=()
migration_script_before=()
skip_tables=false

# PARSE VARIABLE ARGUMENTS WITH getops
options=":c:s:i:u:r:g:tdbkvhao:p:"
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
            migration_script+=("$OPTARG")
            ;;
        p )
            run_script_before=true
            migration_script_before+=("$OPTARG")
            ;;
        o )
            run_script=true
            migration_script+=("$OPTARG")
            ;;
        t )
            tables=true
            ;;
        b )
            tables=false
            skip_tables=true
            ;;
        r )
            ren_app=true
            ;;
		g )
			git_branch=$OPTARG
            ;;
        c )
            conf=$OPTARG
            ;;
        k )
            docker=true
            restart_apache=false
            ;;
        a )
            restart_apache=false
            ;;
        h )
            usage
            exit 1
            ;;
        v )
            echo $APP_VERSION
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

operation="install"
operation_scope="$install_type"
if [ $upgrade == true ]; then
    operation="upgrade"
    operation_scope="$upgrade_type"
fi

# Default to loading initial tables on installs unless explicitly skipped.
if [ "$operation" = "install" ] && [ "$skip_tables" = false ] && [ "$tables" = false ]; then
    tables=true
fi

load_install_config
checkout_git_revision
user=${SUDO_USER:-$USER}
check_requirements

if [[ "$operation_scope" == "full" || "$operation_scope" == "dep" ]]; then
    run_dependency_stage "$operation"
    if [ "$operation_scope" = "dep" ]; then
        log_info "Dependency stage completed."
        exit 0
    fi
fi

if [[ "$operation_scope" == "full" || "$operation_scope" == "app" ]]; then
    if [ "$operation" = "install" ]; then
        install_application_files
    else
        upgrade_application_files
    fi
    if [ $restart_apache == true ]; then
        restart_apache_service
    fi
    exit 0
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
