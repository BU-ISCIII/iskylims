# iSkyLIMS description
[![Django](https://img.shields.io/static/v1?label=Django&message=3.2.8&color=blue?style=plastic&logo=django)](https://github.com/django/django)
[![Python](https://img.shields.io/static/v1?label=Python&message=3.8.10&color=green?style=plastic&logo=Python)](https://www.python.org/)
[![Bootstrap](https://img.shields.io/badge/Bootstrap-v3.0-blueviolet?style=plastic&logo=Bootstrap)](https://getbootstrap.com)
[![version](https://img.shields.io/badge/version-2.2.2-orange?style=plastic&logo=GitHub)](https://github.com/BU-ISCIII/iskylims.git)

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

## Install iSkyLIMS in your server running ubuntu
 
### Pre-requesites
Before starting the installation check :
-   You have **sudo privileges** to install the additional software packets that iSkyLIMS needs.
-   Database (MySQL/MariaDB) is running  
-   Local server configured for sending emails
-   Apache server is running on local server
-   Python installed in server
-   Dependencies:
     - lsb_release:
     RedHat/CentOS: ```yum install redhat-lsb-core```
     Ubuntu: ```xxx`

#### Clone github repository
Open a linux terminal and move to a directory where relecov code will be 
downloaded
```bash
cd <your personal folder>
git clone https://github.com/BU-ISCIII/iSkyLIMS.git iSkyLIMS
cd iSkyLIMS
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

iSkyLIMS should be installed on the "/opt" directory. Before start the installation be sure you have sudo priveleges.

Execute the following commands in a linux terminal.

```bash

sudo bash install.sh
```

After installation is completed open you navigator typing "localhost" or the "server local IP".

# iSkyLIMS documentation

iSkyLIMS documentation is available at https://iskylims.readthedocs.io/en/latest

