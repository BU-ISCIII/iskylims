# iSkyLIMS Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [3.1.0dev] - 2025-XX-XX : <https://github.com/BU-ISCIII/iskylims/releases/tag/3.1.0>

### Credits

- [Sara Monzón](https://github.com/saramonzon)
- [Luis Chapado](https://github.com/luissian)
- [Daniel Valle](https://github.com/Daniel-VM)
- [Pablo Mata](https://github.com/Shettland)
- [Sergio Olmos]()

#### Added Enhancements

- Added setting for HTTPS forwarding [#257](https://github.com/BU-ISCIII/iskylims/pull/257)
- Allow switching to Git SHA or Version Tag and restore initial state [#270](https://github.com/BU-ISCIII/iskylims/pull/270)
- Created graphics for the services that were re-analyzed [#290](https://github.com/BU-ISCIII/iskylims/pull/290)
- Enhance both `install.sh` and `docker_install.sh`, fix data loading issues [#327](https://github.com/BU-ISCIII/iskylims/pull/327)

#### Fixes

- Fixed minor bugs and improved README [#252](https://github.com/BU-ISCIII/iskylims/pull/252)
- Fixed issue #283: Error in Services Statistics per classification area [#286](https://github.com/BU-ISCIII/iskylims/pull/286)
- Fixed issue #289 [#291](https://github.com/BU-ISCIII/iskylims/pull/291)
- Fixed installation script issue where logs symbolic link was not created if it already existed (#256) [#262](https://github.com/BU-ISCIII/iskylims/pull/262)
- Fixed excessive email notifications during crontab process (#266) [#262](https://github.com/BU-ISCIII/iskylims/pull/262)
- Fixed issue where sample names could not be repeated, making SampleID unique (#26) [#262](https://github.com/BU-ISCIII/iskylims/pull/262)
- Prevent underscores in sample names (#73) [#262](https://github.com/BU-ISCIII/iskylims/pull/262)
- Fixed incorrect ordering of service states in `first_install_tables.json` (#265) [#267](https://github.com/BU-ISCIII/iskylims/pull/267)
- Fixed incorrect confirmation email text after resolution (#261) [#267](https://github.com/BU-ISCIII/iskylims/pull/267)
- Fixed issue where users couldn't search service/project by sample name (#264) [#267](https://github.com/BU-ISCIII/iskylims/pull/267)
- Fixed AttributeError when no username is found in wetlab project (#250) [#267](https://github.com/BU-ISCIII/iskylims/pull/267)
- Fixed issue where barcode count conversion to integer failed (#158) [#267](https://github.com/BU-ISCIII/iskylims/pull/267)
- Fixed deletion issue where removing a run also deleted pools and library preparations (#180) [#267](https://github.com/BU-ISCIII/iskylims/pull/267)
- Fixed deprecated `STATUS_CHOICES` usage in Django versions higher than 3.1.x (#263) [#267](https://github.com/BU-ISCIII/iskylims/pull/267)
- Fixed issue where services could not be searched by service type (#78) [#267](https://github.com/BU-ISCIII/iskylims/pull/267)
- Fixed issue [#338](https://github.com/BU-ISCIII/iskylims/issues/338): Removed unnecessary hidden input passing a large JSON object, now using session storage [#344](https://github.com/BU-ISCIII/iskylims/pull/344)
- Fixed email error handling for multiple notification types. [#346](https://github.com/BU-ISCIII/iskylims/pull/346)
- Replaced all references to `Molecule Code ID` with `Extraction Code ID` for consistency.
- Improved message display in Handling Library Preparation.
- Fixed incomplete code execution when storing protocol values. (#349) [#352](https://github.com/BU-ISCIII/iskylims/pull/352)
- Increased maximum length for `prefix_protocol` to prevent data errors. (#350) [#352](https://github.com/BU-ISCIII/iskylims/pull/352)
- Corrected exception handling, replacing incorrect exception type with `AttributeError`. (#351) [#352](https://github.com/BU-ISCIII/iskylims/pull/352)
- Fixed DataError - Value Too Long for prefix_protocol #350: Increased lenght for field prefix_protocol [#352](https://github.com/BU-ISCIII/iskylims/pull/352)

#### Changed

- Updated installation script with variable modules for more flexibility [#269](https://github.com/BU-ISCIII/iskylims/pull/269)
- Updated installation script to remove commas in values of `rawtobunbarcode` table [#277](https://github.com/BU-ISCIII/iskylims/pull/277)
- Updated installation documentation and script, fixing small issues [#284](https://github.com/BU-ISCIII/iskylims/pull/284)
- Unify main and develop branches [#334](https://github.com/BU-ISCIII/iskylims/pull/334)
- Increased max upload memory size. (#328)
- Renamed method `get_delivery_date` to `get_delivered_date` for clarity.
- Improved query performance and excluded rejected/archived services from ongoing list. (#299)
- Updated sample metadata fields with standardized ontology mappings and schema alignment.[#358](https://github.com/BU-ISCIII/iskylims/pull/358)

#### Removed

- Dummy fix in usage line [#271](https://github.com/BU-ISCIII/iskylims/pull/271)

#### Requirements

| Package             | Last release Version  | New release Version  |
|:--------------------|:---------------------|:----------------------------|
| wheel               | 0.37.1               | 0.44.0                      |
| asn1crypto          | 1.5.0                | 1.5.1                       |
| bcrypt              | 4.0.1                | 4.2.0                       |
| biopython           | 1.79                 | 1.84                        |
| cryptography        | 38.0.3               | 43.0.1                      |
| Django              | 4.2                  | 4.2.15                      |
| django-crispy-forms | 2.0                  | 2.3                         |
| django-js-asset     | 2.0.0                | 2.2.0                       |
| django-mptt         | 0.14.0               | 0.16.0                      |
| django-mptt-admin   | 2.4.1                | 2.6.2                       |
| django-cleanup      | 7.0.0                | 8.1.0                       |
| mod_wsgi            | 4.9.4                | 5.0.0                       |
| mysqlclient         | 2.0.3                | 2.2.6                       |
| paramiko            | 3.1.0                | 3.4.1                       |
| jsonschema          | 4.17.3               | 4.23.0                      |
| django_extensions   | 3.2.1                | 3.2.3                       |
| djangorestframework | 3.14.0               | 3.15.2                      |
| drf-yasg            | 1.21.5               | 1.21.7                      |
| pandas              | 1.5.3                | 2.2.2                       |
| openpyxl            | 3.1.1                | 3.1.5                       |
| setuptools          |                      | 75.2.0                      |
