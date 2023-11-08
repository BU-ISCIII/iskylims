# iSkyLIMS

[![Django](https://img.shields.io/static/v1?label=Django&message=4.2&color=blue?style=plastic&logo=django)](https://github.com/django/django)
[![Python](https://img.shields.io/static/v1?label=Python&message=3.8.10&color=green?style=plastic&logo=Python)](https://www.python.org/)
[![Bootstrap](https://img.shields.io/badge/Bootstrap-v5.0-blueviolet?style=plastic&logo=Bootstrap)](https://getbootstrap.com)
[![version](https://img.shields.io/badge/version-3.0.0-orange?style=plastic&logo=GitHub)](https://github.com/BU-ISCIII/iskylims.git)

The introduction of massive sequencing (MS) in genomics facilities has meant an exponential growth in data generation, requiring a precise tracking system, from library preparation to fastq file generation, analysis and delivery to the researcher. Software designed to handle those tasks are called Laboratory Information Management Systems (LIMS), and its software has to be adapted to their own genomics laboratory particular needs. iSkyLIMS is born with the aim of helping with the wet laboratory tasks, and implementing a workflow that guides genomics labs on their activities from library preparation to data production, reducing potential errors associated to high throughput technology, and facilitating the quality control of the sequencing. Also, iSkyLIMS connects the wet lab with dry lab facilitating data analysis by bioinformaticians.

![Image](https://github.com/BU-ISCIII/iskylims/blob/main/img/iskylims_scheme.png)

According to existent infrastructure sequencing is performed on an Illumina NextSeq instrument. Data is stored in NetApp mass storage device and fastq files are generated (bcl2fastq) on a Sun Grid Engine High Performance Computing cluster (SGE-HPC).
Application servers run web applications for bioinformatics analysis (GALAXY), the iSkyLIMS app, and host the MySQL information tier. iSkyLIMS WetLab workflow deals with sequencing run tracking and statistics. Run tracking passes through five states: "recorded” genomics user record the new sequencing run into the system, the process will wait till run is completed by the machine and data is transferred to the mass storage device; “Sample sheet sent” sample sheet file with the sequencing run information will be copied to the run folder for bcl2fastq process; “Processing data” run parameters files are processed and data is stored in the database; “Running stats” demultiplexing data generated in bcl2fastq process is processed and stored into the database, “Completed” all data is processed and stored successfully. Statistics per sample, per project, per run and per investigation are provided, as well as annual and monthly reports. iSkyLIMS DryLab workflow deals with bioinformatics services request and statistics. User request services that can be associated with a sequencing run. Stats and services tracking is provided.

- [iSkyLIMS](#iskylims)
  - [Installation](#installation)
    - [Pre-requisites](#pre-requisites)
    - [iSkyLIMS docker installation](#iskylims-docker-installation)
    - [Install iSkyLIMS in your server running ubuntu/CentOS](#install-iskylims-in-your-server-running-ubuntucentos)
      - [Clone github repository](#clone-github-repository)
      - [Create iskylims database and grant permissions](#create-iskylims-database-and-grant-permissions)
      - [Configuration settings](#configuration-settings)
      - [Run installation script](#run-installation-script)
    - [Upgrade to iSkyLIMS version 3.0.0](#upgrade-to-iskylims-version-300)
      - [Pre-requisites](#pre-requisites-1)
      - [Clone github repository](#clone-github-repository-1)
      - [Configuration settings](#configuration-settings-1)
      - [Running upgrade script](#running-upgrade-script)
        - [Steps requiring root](#steps-requiring-root)
        - [Steps not requiring root](#steps-not-requiring-root)
      - [What to do if something fails](#what-to-do-if-something-fails)
    - [Final configuration steps](#final-configuration-steps)
      - [SAMBA configurarion](#samba-configurarion)
      - [Email verification](#email-verification)
      - [Configure Apache server](#configure-apache-server)
      - [Verification of the installation](#verification-of-the-installation)
    - [iSkyLIMS documentation](#iskylims-documentation)

## Installation

For any problems or bug reporting please post us an [issue](https://github.com/BU-ISCIII/iSkyLIMS/issues)

### Pre-requisites

Before starting the installation make sure :

- You have **sudo privileges** to install the additional software packets that iSkyLIMS needs.
- Database MySQL > 8.0 or MariaDB > 10.4
- Local server configured for sending emails
- Apache server v2.4
- git > 2.34
- Python > 3.8
- Connection to samba shared folder where run folders are stored (p.e galera/NGS_Data)
- Dependencies:
  - lsb_release:
    - RedHat/CentOS: ```yum install redhat-lsb-core```
    - Ubuntu: ```apt install lsb-core lsb-release```

### iSkyLIMS docker installation

You can test iSkyLIMS by creating a docker container on your local machine.

Clone the iSkyLIMS github repository and run the docker script to create the docker

```bash
git clone https://github.com/BU-ISCIII/iSkyLIMS.git iSkyLIMS
sudo bash docker_install.sh
```

The script creates a docker compose container with 3 services:

- web1: contains the iSkyLIMS web application
- db1: contains the mySQL database
- samba: contains samba server

After Docker is created and services are up, database structure and initial data are loaded into database. When this step is completed, you will be asked to define the super user which will have access to django admin pages. You can type any name, but we recommend that you use "admin", because admin user is requested later on when defining the initial settings.

Follow the prompt message to create the super user account.

When script ends open your navigator typing **localhost:8001** to access to iSkyLIMS

### Install iSkyLIMS in your server running ubuntu/CentOS

#### Clone github repository

Open a linux terminal and move to a directory where iSkyLIMS code will be
downloaded

```bash
cd < your personal folder >
git clone https://github.com/BU-ISCIII/iskylims.git iskylims
cd iskylims
```

#### Create iskylims database and grant permissions

1. Create a new database named "iskylims" (this is mandatory)
2. Create a new user with permission to read and modify that database.
3. Write down user, passwd and db server info.

#### Configuration settings

Copy the initial setting template into a file named install_settings.txt

```bash
cp conf/template_install_settings.txt install_settings.txt
```

Open with your favourite editor the configuration file to set your own values for
database ,email settings and the local IP of the server where iSkyLIMS will run.

```bash
nano install_settings.txt
```

#### Run installation script

iSkyLIMS should be installed on the "/opt" directory.

You will need sudo privileges for installing dependencies. In order to handle different installation responsibilities inside the organization, where you may not be the person with root privileges, our instalation script has these options in ```--install``` parameter:

- dep: to install the software packages as well as python packages inside the virtual environment. Root is needed.
- app: to install only the iSkyLIMS application software without need of being root.
- full: if you directly have root permissions you can install both deps and app at the same time with this option.

Execute one of the following commands in a linux terminal to install, according as
above description.

```bash
# to install only software packages dependences
sudo bash install.sh --install dep

# to install only iSkyLIMS application
bash install.sh --install app

# to install both software
sudo bash install.sh --install full
```

### Upgrade to iSkyLIMS version 3.0.0

If you have already iSkyLIMS on version 2.3.0 you can upgrade to the latest stable version 3.0.0.

Version 3.0.0 is a major release with important upgrades in third parties dependencies like bootstrap. Also, we 've done a huge work on refactoring and variables/function renaming that affects the database. For more details about the changes see the release notes.

#### Pre-requisites

Because in this upgrade many tables in database are modified it is required that you backup:

- iSkyLIMS database
- iSkyLIMS folder (complete installation folder, p.e /opt/iSkyLIMS)

It is highly recomended that you made these backups and keep them safely in case of upgrade failure, to recover your system.

#### Clone github repository

We've also change the way that iSkyLIMS is installed and upgraded. From now on iskylims is downloaded in a user folder and installed elsewhere (p.e /opt/).

Open a linux terminal and move to a directory where iSkyLIMS code will be
downloaded

```bash
cd < your personal folder >
git clone https://github.com/BU-ISCIII/iSkyLIMS.git iskylims
cd iskylims
```

#### Configuration settings

Copy the initial setting template into a file named install_settings.txt

```bash
cp conf/template_install_settings.txt install_settings.txt
```

Open with your favourite editor the configuration file to set your own values for
database ,email settings and the local IP of the server where iSkyLIMS will run.
> If you use a windows-based system for modifying the file, make sure the file is saved using a linux-friendly encoding like ASCII or UTF-8

```bash
sudo nano install_settings.txt
```

#### Running upgrade script

If your organization requires that dependencies / stuff that needs root are installed by a different person that install the application the you can use the install script in several steps as follows.

First you need to rename the folder app name in the installation folder (`/opt/iSkyLIMS`):

##### Steps requiring root

```bash
# You need root for this operation
sudo mv /opt/iSkyLIMS /opt/iskylims
```

Make sure that the installation folder has the correct permissions so the person installing the app can write in that folder.

```bash
# In case you have a script for this task. You'll need to adjust this script according to the name changing: /opt/iSkyLIMS to /opt/iskylims
/scripts/hardening.sh
```

In the linux terminal execute one of the following command that fit better to you:

```bash
# to upgrade only software packages dependences. NEEDS ROOT.
sudo bash install.sh --upgrade dep

# to install both software. NEEDS ROOT.
sudo bash install.sh --upgrade full  --ren_app --script drylab_service_state_migration --script rename_app_name --script rename_sample_sheet_folder --script migrate_sample_type --script  migrate_optional_values --tables
```

##### Steps not requiring root

Next you need to upgrade iskylims app. Please use the command below:

```bash
# to upgrade only iSkyLIMS application including changes required in this release. DOES NOT NEED ROOT.
bash install.sh --upgrade app --ren_app --script drylab_service_state_migration --script rename_app_name --script rename_sample_sheet_folder --script migrate_sample_type --script  migrate_optional_values --tables
```

Make sure that the installation folder has the correct permissions.

```bash
# In case you have a script for this task. Some paths have changed in this version, so you may need to adjust your hardening script.
/scripts/hardening.sh
```

#### What to do if something fails

When we upgrade using the installation script we are performing several changes in the database. If something fails we need to restore the app situation before anything happened and start all over.

We need to copy back the full `/opt/iSkyLIMS` folder back to `/opt` (or your installation path preference), and restore the database doing something like this:

```bash
sudo rm -rf /opt/iskylims
sudo cp -r /home/dadmin/backup_prod/iSkyLIMS/ /opt/
sudo /scripts/hardening.sh
mysql -u iskylims -h dmysqlps.isciiides.es
# drop database iskylims;
# create database iskylims;
mysql -u iskylims -h dmysqlps.isciiides.es iskylims < /home/dadmin/backup_prod/bk_iSkyLIMS_202310160737.sql
```

### Final configuration steps

#### SAMBA configurarion

- Login with admin account.
- Go to Massive sequencing
![go_to_wetlab](img/got_to_wetlab.png){width:50px}
- Go to Configuration -> Samba configuration
- Fill the form with the appropiate params for the samba shared folder:
![samba form](img/samba_form.png)

#### Email verification

- Go to Massive sequencing
- Go to Configuration -> Email configuration
- Fill the form with the needed params for your email configuration and try to send a test email.

#### Configure Apache server

Copy the apache configuration file according to your distribution inside the apache configutation directory and rename it to iskylims.conf

#### Verification of the installation

Open the navigator and type "localhost" or the "server local IP" and check that iSkyLIMs is running.

You can also check some of the functionality, while also checking samba and database connections using:

- Go to [configuration test](https://iskylims.isciii.es/wetlab/configurationTest/)
- Click submit
- Check all tabs so every connectin is successful.
- Run the 3 tests for each sequencing machine: MiSeq, NextSeq and NovaSeq.

### iSkyLIMS documentation

iSkyLIMS documentation is available at [https://iskylims.readthedocs.io/en/latest](https://iskylims.readthedocs.io/en/latest)
