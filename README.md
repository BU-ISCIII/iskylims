# iSkyLIMS

[![Django](https://img.shields.io/static/v1?label=Django&message=4.2&color=blue?style=plastic&logo=django)](https://github.com/django/django)
[![Python](https://img.shields.io/static/v1?label=Python&message=3.8.10&color=green?style=plastic&logo=Python)](https://www.python.org/)
[![Bootstrap](https://img.shields.io/badge/Bootstrap-v5.0-blueviolet?style=plastic&logo=Bootstrap)](https://getbootstrap.com)
[![version](https://img.shields.io/badge/version-3.0.0-orange?style=plastic&logo=GitHub)](https://github.com/BU-ISCIII/iskylims.git)

The introduction of massive sequencing (MS) in genomics facilities has meant an exponential growth in data generation, requiring a precise tracking system, from library preparation to fastq file generation, analysis and delivery to the researcher. Software designed to handle those tasks are called Laboratory Information Management Systems (LIMS), and its software has to be adapted to their own genomics laboratory particular needs. iSkyLIMS is born with the aim of helping with the wet laboratory tasks, and implementing a workflow that guides genomics labs on their activities from library preparation to data production, reducing potential errors associated to high throughput technology, and facilitating the quality control of the sequencing. Also, iSkyLIMS connects the wet lab with dry lab facilitating data analysis by bioinformaticians.

![Image](img/iskylims_scheme.png)

According to existent infrastructure sequencing is performed on an Illumina NextSeq instrument. Data is stored in NetApp mass storage device and fastq files are generated (bcl2fastq) on a Sun Grid Engine High Performance Computing cluster (SGE-HPC).
Application servers run web applications for bioinformatics analysis (GALAXY), the iSkyLIMS app, and host the MySQL information tier. iSkyLIMS WetLab workflow deals with sequencing run tracking and statistics. Run tracking passes through five states: "recorded” genomics user record the new sequencing run into the system, the process will wait till run is completed by the machine and data is transferred to the mass storage device; “Sample sheet sent” sample sheet file with the sequencing run information will be copied to the run folder for bcl2fastq process; “Processing data” run parameters files are processed and data is stored in the database; “Running stats” demultiplexing data generated in bcl2fastq process is processed and stored into the database, “Completed” all data is processed and stored successfully. Statistics per sample, per project, per run and per investigation are provided, as well as annual and monthly reports. iSkyLIMS DryLab workflow deals with bioinformatics services request and statistics. User request services that can be associated with a sequencing run. Stats and services tracking is provided.

- [iSkyLIMS](#iskylims)
  - [Get the code (required)](#get-the-code-required)
  - [Choose your path](#choose-your-path)
  - [Docker deployment](#docker-deployment)
    - [Local test stack](#local-test-stack)
    - [Production container](#production-container)
    - [Upgrade docker deployment](#upgrade-docker-deployment)
  - [Bare-metal deployment (Ubuntu/CentOS)](#bare-metal-deployment-ubuntucentos)
    - [Install](#install)
      - [Prerequisites](#prerequisites)
      - [Clone the repository](#clone-the-repository)
      - [Prepare the database](#prepare-the-database)
      - [Configure install\_settings.txt](#configure-install_settingstxt)
      - [Run install.sh](#run-installsh)
    - [Upgrade (3.0.x to 3.1.x)](#upgrade-30x-to-31x)
      - [Back up first](#back-up-first)
      - [Refresh code and settings](#refresh-code-and-settings)
      - [Run upgrade steps requiring root](#run-upgrade-steps-requiring-root)
      - [Run upgrade steps without root](#run-upgrade-steps-without-root)
  - [What to do if something fails](#what-to-do-if-something-fails)
  - [Final configuration steps](#final-configuration-steps)
    - [SAMBA configurarion](#samba-configurarion)
    - [Email verification](#email-verification)
    - [Configure Apache server](#configure-apache-server)
    - [Verification of the installation](#verification-of-the-installation)
  - [iSkyLIMS documentation](#iskylims-documentation)

For any problems or bug reporting please post us an [issue](https://github.com/BU-ISCIII/iSkyLIMS/issues)

## Get the code (required)

All installation paths assume you already cloned the repository:

```bash
git clone https://github.com/BU-ISCIII/iskylims.git iskylims
cd iskylims
```

## Choose your path

- **Docker (local test)**: spin up MySQL + Samba + iSkyLIMS with demo data to try the app quickly.
- **Docker (production container)**: deploy only the application container, pointing to your existing DB/Samba.
- **Bare-metal**: install or upgrade directly on Ubuntu/CentOS hosts with `install.sh`.

## Docker deployment

Prerequisites for Docker-based installs:

- Docker Engine + Docker Compose v2
- git (to clone the repository)

### Local test stack

Bring up a full test stack (database, Samba, app) plus fixtures and demo data:

```bash
bash docker_install.sh --test
```

This uses `docker-compose.test.yml` by default.

Defaults can be customised:

- `--demo_data /path/to/iskylims_demo_data.tar.gz` to reuse a local demo archive (otherwise it is downloaded).
- `--skip_demo_data` or `--skip_test_data` to avoid loading extra data.
- `--install_type` (`full` by default) and `--git_revision` to control the build.
- `--script` to run one or more Django migration scripts through `install.sh` (repeat the flag as needed).

Example running a migration script during Docker install:

```bash
bash docker_install.sh --test --script migrate_optional_values
```

When the script finishes, open `http://localhost:8001` and follow the prompt to create the Django superuser.

### Production container

Deploy the iSkyLIMS container against external MySQL/Samba services:

1. Copy and edit the production settings template:

    ```bash
    cp conf/docker_production_settings.txt conf/my_prod_settings.txt
    # edit conf/my_prod_settings.txt with your DB/Samba details
    ```

2. Build and run in production mode (uses `docker-compose.prod.yml` by default):

    ```bash
    bash docker_install.sh --install_conf conf/my_prod_settings.txt
    ```

   Use `--compose_file` to override the compose file or `--install_type`/`--git_revision` to change the build.

3. If this is a fresh install, create the Django superuser when prompted and complete the Samba configuration in the UI.

### Upgrade docker deployment

Re-deploy the application container against an existing production database without touching data:

```bash
bash docker_install.sh --install_conf conf/my_prod_settings.txt --action upgrade
```

The upgrade path rebuilds/restarts the container and runs `install.sh` inside the app container, which regenerates migrations, applies them with `--fake-initial`, and skips superuser/demo/test data loading.

## Bare-metal deployment (Ubuntu/CentOS)

### Install

#### Prerequisites

- **sudo privileges** for dependency installation
- MySQL > 8.0 or MariaDB > 10.4
- Apache 2.4
- git > 2.34
- Python > 3.8
- Local email sender configured
- Access to the Samba share where run folders live
- `lsb_release` package (`yum install redhat-lsb-core` on RedHat/CentOS, `apt install lsb-core lsb-release` on Ubuntu)

#### Clone the repository

```bash
cd <your working directory>
git clone https://github.com/BU-ISCIII/iskylims.git iskylims
cd iskylims
```

#### Prepare the database

1. Create a database named `iskylims`.
2. Create a user with read/write permissions on that database.
3. Note the database host, port, user, and password for the settings file.

#### Configure install_settings.txt

```bash
cp conf/template_install_settings.txt install_settings.txt
nano install_settings.txt
```

Set your database, email, server IP/URL, and logging preferences in that file.

#### Run install.sh

iSkyLIMS is installed to `/opt/iskylims` by default. The single `install.sh` script handles both dependencies and the app; choose what you need with `--install`:

- `dep`: install system and Python dependencies (requires sudo).
- `app`: deploy iSkyLIMS code, update settings, run migrations, and collect static files (no sudo needed).
- `full`: run both stages in sequence.

Examples:

```bash
# only software dependencies
sudo bash install.sh --install dep

# only iSkyLIMS application
bash install.sh --install app --git_revision main --tables

# dependencies + application
sudo bash install.sh --install full --git_revision main --tables
```

- Add `--tables` to load the initial fixtures on first-time installs, or `--skip_tables` if you want to skip them.
- Capture logs for troubleshooting with `tee`:

  ```bash
  sudo bash install.sh --install full --git_revision main --tables 2>&1 | tee install_full.log
  ```

- If Apache is managed elsewhere, skip the automatic restart with `--skip_apache_restart`.

### Upgrade (3.0.x to 3.1.x)

Follow these steps to move from version 3.0.0 to the 3.1.x series.

#### Back up first

- Full backup of the `iskylims` database.
- Full backup of the installation folder (for example `/opt/iskylims`).
- If you use library pools, export them before upgrading:

  ```bash
  mysql --user=<db_user> --password=<db_password> --host=<db_server_ip> --port=<db_port> iskylims \
    -e "SELECT * FROM wetlab_library_pool" > <backup_folder>/backup_lib_pool.sql
  ```

#### Refresh code and settings

```bash
cd <your working directory>/iskylims
git pull
cp conf/template_install_settings.txt install_settings.txt
sudo nano install_settings.txt
```

Ensure the file uses Linux-friendly encoding (UTF-8/ASCII) if you edit it on Windows.

#### Run upgrade steps requiring root

Update system and Python dependencies:

```bash
sudo bash install.sh --upgrade dep 2>&1 | tee install_full.log
```

Make sure the installation directory permissions allow the non-root step to write to `/opt/iskylims` (adapt your hardening script if paths changed).

#### Run upgrade steps without root

Upgrade the application code and database:

```bash
# with library pool restore
bash install.sh --upgrade app --script <backup_folder>/backup_lib_pool.sql --git_revision main --tables

# without library pool restore
bash install.sh --upgrade app --git_revision main --tables

# example running a migration script during upgrade
bash install.sh --upgrade app --script migrate_optional_values --git_revision main --tables
```

Or run everything in one go:

```bash
sudo bash install.sh --upgrade full --git_revision main --tables
```

Upgrades regenerate migrations and apply them with `--fake-initial` so existing tables remain intact, matching the Docker workflow.

## What to do if something fails

When we upgrade using the installation script we are performing several changes in the database. If something fails we need to restore the app situation before anything happened and start all over.

We need to copy back the full `/opt/iskylims` folder back to `/opt/iskylims` (or your installation path preference), and restore the database doing something like this:

```bash
sudo rm -rf /opt/iskylims
sudo cp -r /home/dadmin/backup_prod/iSkyLIMS/ /opt/
sudo /scripts/hardening.sh
mysql -u iskylims -p -h dmysqlps.isciiides.es
# drop database iskylims;
# create database iskylims;
mysql -u iskylims -p -h dmysqlps.isciiides.es iskylims < /home/dadmin/backup_prod/bk_iSkyLIMS_202310160737.sql
```

## Final configuration steps

### SAMBA configurarion

- Login with admin account.
- Go to Massive sequencing
![go_to_wetlab](img/got_to_wetlab.png){width:50px}
- Go to Configuration -> Samba configuration
- Fill the form with the appropiate params for the samba shared folder:
![samba form](img/samba_form.png)

### Email verification

- Go to Massive sequencing
- Go to Configuration -> Email configuration
- Fill the form with the needed params for your email configuration and try to send a test email.

### Configure Apache server

Copy the apache configuration file according to your distribution inside the apache configutation directory and rename it to iskylims.conf

### Verification of the installation

Open the navigator and type "localhost" or the "server local IP" and check that iSkyLIMs is running.

You can also check some of the functionality, while also checking samba and database connections using:

- Go to [configuration test](https://iskylims.isciii.es/wetlab/configurationTest/)
- Click submit
- Check all tabs so every connectin is successful.
- Run the 3 tests for each sequencing machine: MiSeq, NextSeq and NovaSeq.

## iSkyLIMS documentation

iSkyLIMS documentation is available at [https://iskylims.readthedocs.io/en/latest](https://iskylims.readthedocs.io/en/latest)
