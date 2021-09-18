# iSkyLIMS description

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

Then open your navigator and type, localhost:8000

## Install iSkyLIMS in your server running ubuntu
### Pre-requesites

You need to be sudo priviledges to install 
#### Database

#### email

#### Configuration settings

### Run installation script


```bash
git clone https://github.com/BU-ISCIII/iSkyLIMS.git iSkyLIMS
bash install.sh
```

# iSkyLIMS documentation
iSkyLIMS documentation is available at https://iskylims.readthedocs.io/en/latest

