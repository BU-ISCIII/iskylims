# iSkyLIMS

[![Django](https://img.shields.io/static/v1?label=Django&message=3.2.8&color=blue?style=plastic&logo=django)](https://github.com/django/django)
[![Python](https://img.shields.io/static/v1?label=Python&message=3.8.10&color=green?style=plastic&logo=Python)](https://www.python.org/)
[![Bootstrap](https://img.shields.io/badge/Bootstrap-v5.0-blueviolet?style=plastic&logo=Bootstrap)](https://getbootstrap.com)
[![version](https://img.shields.io/badge/version-2.3.1-orange?style=plastic&logo=GitHub)](https://github.com/BU-ISCIII/iskylims.git)

The introduction of massive sequencing (MS) in genomics facilities has meant an exponential growth in data generation, requiring a precise tracking system, from library preparation to fastq file generation, analysis and delivery to the researcher. Software designed to handle those tasks are called Laboratory Information Management Systems (LIMS), and its software has to be adapted to their own genomics laboratory particular needs. iSkyLIMS is born with the aim of helping with the wet laboratory tasks, and implementing a workflow that guides genomics labs on their activities from library preparation to data production, reducing potential errors associated to high throughput technology, and facilitating the quality control of the sequencing. Also, iSkyLIMS connects the wet lab with dry lab facilitating data analysis by bioinformaticians.

<img src="https://github.com/BU-ISCIII/iSkyLIMS/blob/main/img/iskylims_scheme.png" width="900">

According to existent infrastructure sequencing is performed on an Illumina NextSeq instrument. Data is stored in NetApp mass storage device and fastq files are generated (bcl2fastq) on a Sun Grid Engine High Performance Computing cluster (SGE-HPC).
Application servers run web applications for bioinformatics analysis (GALAXY), the iSkyLIMS app, and host the MySQL information tier. iSkyLIMS WetLab workflow deals with sequencing run tracking and statistics. Run tracking passes through five states: "recorded” genomics user record the new sequencing run into the system, the process will wait till run is completed by the machine and data is transferred to the mass storage device; “Sample sheet sent” sample sheet file with the sequencing run information will be copied to the run folder for bcl2fastq process; “Processing data” run parameters files are processed and data is stored in the database; “Running stats” demultiplexing data generated in bcl2fastq process is processed and stored into the database, “Completed” all data is processed and stored successfully. Statistics per sample, per project, per run and per investigation are provided, as well as annual and monthly reports. iSkyLIMS DryLab workflow deals with bioinformatics services request and statistics. User request services that can be associated with a sequencing run. Stats and services tracking is provided.

- [iSkyLIMS](#iskylims)
  - [Installation](#installation)
  - [iSkyLIMS docker installation](#iskylims-docker-installation)
  - [Install iSkyLIMS in your server running ubuntu/CentOS](#install-iskylims-in-your-server-running-ubuntucentos)
    - [Pre-requisites](#pre-requisites)
    - [Clone github repository](#clone-github-repository)
    - [Create relecov database and grant permissions](#create-relecov-database-and-grant-permissions)
    - [Configuration settings](#configuration-settings)
    - [Run installation script](#run-installation-script)
    - [Configure Apache server](#configure-apache-server)
    - [Installation verification](#installation-verification)
  - [Upgrade to iSkyLIMS version 3.0.0](#upgrade-to-iskylims-version-300)
    - [Pre-requisites](#pre-requisites-1)
    - [Executing the upgrade](#executing-the-upgrade)
    - [Clone github repository](#clone-github-repository-1)
    - [Configuration settings](#configuration-settings-1)
    - [Runing upgrade script](#runing-upgrade-script)
    - [Verification installation](#verification-installation)
  - [iSkyLIMS documentation](#iskylims-documentation)

## Installation

For any problems or bug reporting please post us an [issue](https://github.com/BU-ISCIII/iSkyLIMS/issues)

## iSkyLIMS docker installation

You can test iSkyLIMS by creating a docker container on your local machine.

Clone the iSkyLIMS github repository and run the docker script to create the docker

```bash
git clone https://github.com/BU-ISCIII/iSkyLIMS.git iSkyLIMS
bash docker_iskylims_install.sh
```

The script creates a docker compose container with 2 services:

- web1: contains the iSkyLIMS web application
- db1: contains the mySQL database

After Docker is created and services are up, database structure and initial data are loaded into database. When this step is completed, you will be asked to define the super user which will have access to django admin pages. You can type any name, but we recommend that you use "admin", because admin user is requested later on when defining the initial settings.

Follow the prompt message to create the super user account.

When script ends open your navigator typing **localhost:8000** to access to iSkyLIMS

## Install iSkyLIMS in your server running ubuntu/CentOS

### Pre-requisites

Before starting the installation make sure :

- You have **sudo privileges** to install the additional software packets that iSkyLIMS needs.
- Database (MySQL/MariaDB) is running.
- Local server configured for sending emails
- Apache server is running.
- Python > 3.8 installed
- Connection to samba dir where run folders are stored.
*- Dependencies:
  - lsb_release:
    - RedHat/CentOS: ```yum install redhat-lsb-core```
    - Ubuntu: ```apt install lsb-core lsb-release```

### Clone github repository

Open a linux terminal and move to a directory where iSkyLIMS code will be
downloaded

```bash
cd <your personal folder>
git clone https://github.com/BU-ISCIII/iskylims.git iskylims
cd iskylims
```

### Create relecov database and grant permissions

1. Create a new database named "iskylims" (this is mandatory)
2. Create a new user with permission to read and modify that database.
3. Write down user, passwd and db server info.

### Configuration settings

Copy the initial setting template into a file named install_settings.txt

```bash
cp conf/template_install_settings.txt install_settings.txt
```

Open with your favourite editor the configuration file to set your own values for
database ,email settings and the local IP of the server where iSkyLIMS will run.

```bash

sudo nano install_settings.txt
```

### Run installation script

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

### Configure Apache server

Copy the apache configuration file according to your distribution inside the apache configutation directory and rename it to iskylims.conf

### Installation verification

After installation is completed and apache server is up and running open you navigator typing "localhost" or the "server local IP".

## Upgrade to iSkyLIMS version 3.0.0

If you have already iSkyLIMS on version 2.3.0 you can upgrade to the latest stable version 3.0.0.

Version 3.0.0 is a major release with important upgrades in third parties dependencies like bootstrap. Also, we 've done a huge work on refactoring and variables/function renaming that affects the database. For more details about the changes see the release notes.

### Pre-requisites

Because in this upgrade many tables in database are modified it is required that you backup:

- iSkyLIMS database
- iSkyLIMS document folder

It is highly recomended that you made these backups and keep them safely in case of
upgrade failure, to recover your system.

### Executing the upgrade

We've also change the way that iSkyLIMS is installed and upgraded. From now on iskylims is downloaded in a user folder and installed elsewhere (p.e /opt/).

### Clone github repository

Open a linux terminal and move to a directory where iSkyLIMS code will be
downloaded

```bash
cd <your personal folder>
git clone https://github.com/BU-ISCIII/iSkyLIMS.git iskylims
cd iskylims
```

### Configuration settings

Copy the initial setting template into a file named install_settings.txt

```bash
cp conf/template_install_settings.txt install_settings.txt
```

Open with your favourite editor the configuration file to set your own values for
database ,email settings and the local IP of the server where iSkyLIMS will run.

```bash

sudo nano install_settings.txt
```

### Runing upgrade script

If your organization requires that dependencies / stuff that needs root are installed by a different person that install the application the you can use the install script in several steps as follows.

First you need to rename the folder app name in the installation folder (`/opt/iSkyLIMS`):

```bash
# You need root for this operation
mv /opt/iSkyLIMS /opt/iskylims
```

Make sure that the installation folder has the correct permissions so the person installing the app can write in that folder.

In the linux terminal execute one of the following command that fit better to you:

```bash
# to upgrade only software packages dependences. NEEDS ROOT.
sudo bash install.sh --upgrade dep

# to upgrade only iSkyLIMS application including changes required in this release. DOES NOT NEED ROOT.
bash install.sh --upgrade app --ren_app --script drylab_service_state_migration --tables

# to install both software. NEEDS ROOT.
sudo bash install.sh --upgrade full  --ren_app --script drylab_service_state_migration --tables
```

### Verification installation

Open the navigator and type "localhost" or the "server local IP" and check that iSkyLIMs is running.

You can also check some of the functionality, while also checking samba and database connections using:

- Go to [configuration test](https://iskylims.isciii.es/wetlab/configurationTest/)
- Click submit
- Check all tabs so every connectin is successful.
- Run the 3 tests for each sequencing machine: MiSeq, NextSeq and NovaSeq.

## iSkyLIMS documentation

iSkyLIMS documentation is available at <https://iskylims.readthedocs.io/en/latest>
