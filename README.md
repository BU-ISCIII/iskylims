# iSkyLIMS description
[![Django](https://img.shields.io/static/v1?label=Django&message=3.2.8&color=blue?style=plastic&logo=django)](https://github.com/django/django)
[![Python](https://img.shields.io/static/v1?label=Python&message=3.8.10&color=green?style=plastic&logo=Python)](https://www.python.org/)
[![Bootstrap](https://img.shields.io/badge/Bootstrap-v5.0-blueviolet?style=plastic&logo=Bootstrap)](https://getbootstrap.com)
[![version](https://img.shields.io/badge/version-2.3.1-orange?style=plastic&logo=GitHub)](https://github.com/BU-ISCIII/iskylims.git)

The introduction of massive sequencing (MS) in genomics facilities has meant an exponential growth in data generation, requiring a precise tracking system, from library preparation to fastq file generation, analysis and delivery to the researcher. Software designed to handle those tasks are called Laboratory Information Management Systems (LIMS), and its software has to be adapted to their own genomics laboratory particular needs. iSkyLIMS is born with the aim of helping with the wet laboratory tasks, and implementing a workflow that guides genomics labs on their activities from library preparation to data production, reducing potential errors associated to high throughput technology, and facilitating the quality control of the sequencing. Also, iSkyLIMS connects the wet lab with dry lab facilitating data analysis by bioinformaticians.

<img src="https://github.com/BU-ISCIII/iSkyLIMS/blob/master/img/iSkyLIMS_scheme.png" width="900">

According to existent infrastructure sequencing is performed on an Illumina NextSeq instrument. Data is stored in NetApp mass storage device and fastq files are generated (bcl2fastq) on a Sun Grid Engine High Performance Computing cluster (SGE-HPC).
Application servers run web applications for bioinformatics analysis (GALAXY), the iSkyLIMS app, and host the MySQL information tier. iSkyLIMS WetLab workflow deals with sequencing run tracking and statistics. Run tracking passes through five states: "recorded” genomics user record the new sequencing run into the system, the process will wait till run is completed by the machine and data is transferred to the mass storage device; “Sample sheet sent” sample sheet file with the sequencing run information will be copied to the run folder for bcl2fastq process; “Processing data” run parameters files are processed and data is stored in the database; “Running stats” demultiplexing data generated in bcl2fastq process is processed and stored into the database, “Completed” all data is processed and stored successfully. Statistics per sample, per project, per run and per investigation are provided, as well as annual and monthly reports. iSkyLIMS DryLab workflow deals with bioinformatics services request and statistics. User request services that can be associated with a sequencing run. Stats and services tracking is provided.

# Installation



For any problems or bug reporting please post us an [issue](https://github.com/BU-ISCIII/iSkyLIMS/issues)

## iSkyLIMS docker installation
You can test iSkyLIMS by creating a docker container on your local machine.

Clone the iSkyLIMS github repository and run the docker script to create the docker


```bash
git clone https://github.com/BU-ISCIII/iSkyLIMS.git iSkyLIMS
bash docker_iskylims_install.sh
```
The script creates a Docker compose container with 2 services:

* web1. Which contains the iSkyLIMS web application
* db1. With the mySQL database

After Docker is created and services are up, database structure and initial data are loaded into database. When this step is completed you will
prompt to define the super user which will be the one to connect to "django admin pages". You can type any name, but we recommend that you use "admin" ,
because admin user is requested later on when defining the initial settings.

Follow the prompt message to create the super user account.

When script ends open your navigator typing **localhost:8000** to access to iSkyLIMS

## Install iSkyLIMS in your server running ubuntu/CentOS
 
### Pre-requesites
Before starting the installation check :
-   You have **sudo privileges** to install the additional software packets that iSkyLIMS needs.
-   Database (MySQL/MariaDB) is running  
-   Local server configured for sending emails
-   Apache server is running on local server
-   Python installed in server
-   Dependencies:
     - lsb_release:
        - RedHat/CentOS: ```yum install redhat-lsb-core```
        - Ubuntu: ```apt install lsb_release```

#### Clone github repository
Open a linux terminal and move to a directory where iSkyLIMS code will be 
downloaded
```bash
cd <your personal folder>
git clone https://github.com/BU-ISCIII/iSkyLIMS.git iskylims
cd iskylims
``` 

#### Create relecov database and grant permissions

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

sudo nano install_settings.txt
```

### Run installation script

iSkyLIMS should be installed on the "/opt" directory. 


Before start the installation be sure you have sudo priveleges.


To handle different installation responsibilities inside the organization, where 
software packets are responsability on different unit that application software.


To support this kind of scenario, instalation script has these options in ```--install```
 parameter.

- dep. To install the software packages as well as python packages inside the virtual environment.
- app. To install only the iSkyLIMS application software.
- full. For those organizations that not require separate installation it wil install both, 
dependecies and iSkyLIMS application.

Execute on of the following commands in a linux terminal to install, according as 
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

Copy the apache configuration file according to your distribution inside the 
apache configutation directory and rename it to iskylims.conf

### Installation verification

After installation is completed and apache server is up and running open you navigator
 typing "localhost" or the "server local IP".

# Upgrade to iSkyLIMS version 2.3.1

If you have already iSkyLIMS on version 2.3.0 you can upgrade to the latest stable 
version 2.3.1.

There are many changes in version 2.3.1 we have upgraded to bootstrap 5, renaming 
tables, update code for a better readability and support, fixing bugs. For a more 
deatils of the changes see the release note for this release.

As previous mention these changes in the code structure of iSkyLIMS, were necessary
 to allow us to incorporate iSkyLIMS with new functionalities and make it easy to
 maintain.


## Requirements

Because in this upgrade many tables in database are impacted it is required that
you backup:

- iskylims database
- iSkyLIMS document folder

It is highly recomended that you made these backups and keep them safely in case of 
upgrade failure, to recover your system.

## Executing the upgrade

We have also change the way that iSkyLIMS is upgraded. In previous releases upgrade
changes where done directly on the directory which iSkyLIMS is running. From these 
release git repository is in the local user folder.


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

if your organization requires that software packages are installed by a different
person that install the software application the upgrade script is also updated
to support this situation.

In the linux terminal execute one of the following command that fit better to you

```bash
# to upgrade only software packages dependences
sudo bash install.sh --upgrade dep

# to upgrade only iSkyLIMS application including changes required in this release
bash install.sh --upgrade app --ren_app --script drylab_service_state_migration --tables

# to install both software
sudo bash install.sh --upgrade full  --ren_app --script drylab_service_state_migration --tables
```

### Verification installation.

Open the navigator and type "localhost" or the "server local IP" and check that
 iSkyLIMs is running.


# iSkyLIMS documentation

iSkyLIMS documentation is available at https://iskylims.readthedocs.io/en/latest

