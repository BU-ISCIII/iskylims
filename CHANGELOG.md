# iSkyLIMS Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [3.1.1dev] - 2026-06-18 : <https://github.com/BU-ISCIII/iskylims/compare/3.1.0...develop>

### Credits

- [Sara Monzón](https://github.com/saramonzon)
- [Enrique Sapena Ventura](https://github.com/ESapenaVentura)
- [Sarai Varona](https://github.com/svarona)

#### Added Enhancements

- Added configurable Apache `/server-status` support, restricted to localhost by default.
- Added generation of supercronic jobs directly from Django `CRONJOBS` settings.
- Added supercronic-compatible cron status and fallback handling for wet-lab scheduled processes.
- Added the `sample-project-values-bulk` API endpoint to retrieve selected project values for multiple RELECOV samples. [#401](https://github.com/BU-ISCIII/iskylims/pull/401)

#### Fixes

- Fixed permission repair so cron directories are also prepared with the expected ownership and mode.
- Fixed container cron status checks when the system `crontab` command is unavailable.
- Fixed Django settings bind path naming in production configuration templates.
- Fixed SampleSheet username resolution when dash-delimited user IDs must match database usernames containing dots. [#396](https://github.com/BU-ISCIII/iskylims/pull/396)
- Fixed website file downloads in production by mounting the shared documents volume read-only in the Apache container. [3bf30dfd](https://github.com/BU-ISCIII/iskylims/commit/3bf30dfd9f65f2760c07f8a2a6da9eded41200db)
- Fixed Apache configuration loading order by assigning explicit ordering prefixes to the generated configuration files.

#### Changed

- Updated production configuration templates with clearer host-path and server-status settings.
- Changed production forwarded-header defaults to HTTPS on port 443.
- Separated Apache and application log paths, including the dedicated RELECOV iSkyLIMS application log path.
- Added local production settings files to `.gitignore`.

#### Removed

#### Requirements

| Package      | Last release Version | New release Version |
|:-------------|:---------------------|:--------------------|
| biopython    | 1.84                 | 1.87                |
| cryptography | 44.0.3               | 48.0.1              |
| Django       | 4.2.28               | 4.2.30              |
| paramiko     | 3.4.1                | 5.0.0               |

## [3.1.0] - 2026-05-27 : <https://github.com/BU-ISCIII/iskylims/releases/tag/3.1.0>

### Credits

- [Sara Monzón](https://github.com/saramonzon)
- [Luis Chapado](https://github.com/luissian)
- [Daniel Valle](https://github.com/Daniel-VM)
- [Pablo Mata](https://github.com/Shettland)
- [Sergio Olmos](https://github.com/OPSergio)

#### Added Enhancements

- Added setting for HTTPS forwarding [#257](https://github.com/BU-ISCIII/iskylims/pull/257)
- Allow switching to Git SHA or Version Tag and restore initial state [#270](https://github.com/BU-ISCIII/iskylims/pull/270)
- Created graphics for the services that were re-analyzed [#290](https://github.com/BU-ISCIII/iskylims/pull/290)
- Enhance both `install.sh` and `docker_install.sh`, fix data loading issues [#327](https://github.com/BU-ISCIII/iskylims/pull/327)
- Included thorough description for wetlab API's update_lab() method [#361](https://github.com/BU-ISCIII/iskylims/pull/361)
- Included new API function lab-request-mapping to get LabRequest fields ontology map [#377](https://github.com/BU-ISCIII/iskylims/pull/377)
- Improved responses in API create-sample-data by adding ERROR messages and data [#377](https://github.com/BU-ISCIII/iskylims/pull/377)
- Added support for `script-before` and `script-after` hooks in install script. [#389](https://github.com/BU-ISCIII/iskylims/pull/389)
- Committed baseline migration files and added migrations for develop changes [#389](https://github.com/BU-ISCIII/iskylims/pull/389)
- Added developer notes on how to create and manage migration files [#389](https://github.com/BU-ISCIII/iskylims/pull/389)
- Added documentation describing migration scripts and their related versions [#389](https://github.com/BU-ISCIII/iskylims/pull/389)
- Enabled Docker internal networking for local test installation [#389](https://github.com/BU-ISCIII/iskylims/pull/389)
- Opened Docker network to allow localhost MySQL connection when required [#389](https://github.com/BU-ISCIII/iskylims/pull/389)
- Added test fixtures and installation updates for the new sequencer and SampleSheet v2 support, including admin test-group assignment and new sequencer bootstrap data. Closes [#388](https://github.com/BU-ISCIII/iskylims/issues/388) [#392](https://github.com/BU-ISCIII/iskylims/pull/392)
- Refactored container installation flow to separate staged application install from runtime bootstrap tasks [#393](https://github.com/BU-ISCIII/iskylims/pull/393)
- Updated Docker image build to stage application files at build time and run bootstrap tasks on container start/upgrade [#393](https://github.com/BU-ISCIII/iskylims/pull/393)
- Replaced container cron runtime with supercronic and improved multi-container install configuration [#391](https://github.com/BU-ISCIII/iskylims/pull/391)
- Added production container support for bind-mounted Django settings and generated Apache configuration files.
- Added production compose environment file generation during container installation.
- Added configurable ServerName and Apache log naming from installation settings.
- Added dedicated Docker network configuration to the production compose file.
- Added configuration examples for container bind mounts in installation settings templates.
- Added explicit preparation and ownership handling for production container bind mounts and named volumes. [2ee02a02](https://github.com/BU-ISCIII/iskylims/commit/2ee02a02c7948702e4f771e9ed6b72598880449e)
- Added the `--action fix-permissions` container installer option to repair host bind mounts and mounted application volumes without rebuilding images or running migrations. [b27f9acc](https://github.com/BU-ISCIII/iskylims/commit/b27f9acc6f6a21f04ed8bc353162f4a9d09411e9)
- Added Apache access to the iSkyLIMS documents volume.

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
- Fixed email error manage for multiple notification types. [#346](https://github.com/BU-ISCIII/iskylims/pull/346)
- Replaced all references to `Molecule Code ID` with `Extraction Code ID` for consistency.
- Improved message display in Manage Library Preparation.
- Fixed incomplete code execution when storing protocol values. (#349) [#352](https://github.com/BU-ISCIII/iskylims/pull/352)
- Increased maximum length for `prefix_protocol` to prevent data errors. (#350) [#352](https://github.com/BU-ISCIII/iskylims/pull/352)
- Corrected exception manage, replacing incorrect exception type with `AttributeError`. (#351) [#352](https://github.com/BU-ISCIII/iskylims/pull/352)
- Fixed DataError - Value Too Long for prefix_protocol #350: Increased lenght for field prefix_protocol [#352](https://github.com/BU-ISCIII/iskylims/pull/352)
- Removed --no-cache from docker_install.sh as it only worked with deprecated docker-compose [#356](https://github.com/BU-ISCIII/iskylims/pull/356)
- Removed unused field sample_project_searchable that was leading to errors during migration [#356](https://github.com/BU-ISCIII/iskylims/pull/356)
- Fixed small spacing errors in docker-compose.yml [#356](https://github.com/BU-ISCIII/iskylims/pull/356)
- Fixed KeyError in project schema loading by using default value for missing 'Downloadable' field.[#358](https://github.com/BU-ISCIII/iskylims/pull/358)
- Wetlab api create_sample_data() also creates lab based on submitting_institution data if present [#360](https://github.com/BU-ISCIII/iskylims/pull/360)
- Wetlab api update_lab() now creates new lab if 'create_if_missing' in request.data [#360](https://github.com/BU-ISCIII/iskylims/pull/360)
- Fixed wetlab API's labrequest.serializer update method to work properly [#361](https://github.com/BU-ISCIII/iskylims/pull/361)
- Adapted update_lab serializer call to new serializer update method [#361](https://github.com/BU-ISCIII/iskylims/pull/361)
- Leave missing submitting_fields as empty string instead of crashing in wetlab.api.create_sample_data [#363](https://github.com/BU-ISCIII/iskylims/pull/363)
- Fixed wetlab API create-sample-data error when submitting_institution fields were not provided [#377](https://github.com/BU-ISCIII/iskylims/pull/377)
- Fixed incorrect git revision propagation to Dockerfile during build [#389](https://github.com/BU-ISCIII/iskylims/pull/389)
- Fixed configuration file accessibility from within Docker container [#389](https://github.com/BU-ISCIII/iskylims/pull/389)
- Fixed disk utilization check to correctly resolve application folder path [#389](https://github.com/BU-ISCIII/iskylims/pull/389)
- Fixed incorrect application folder path resolution in crontab scripts [#389](https://github.com/BU-ISCIII/iskylims/pull/389)
- Fixed samplesheet parsing error [#389](https://github.com/BU-ISCIII/iskylims/pull/389)
- Fixed logger inconsistency that prevented exceptions and error messages from being written to the update crontab log. Closes [#390](https://github.com/BU-ISCIII/iskylims/issues/390) [#392](https://github.com/BU-ISCIII/iskylims/pull/392)
- Fixed wetlab crontab processing for the new sequencer, including run discovery, completion checks, RunInfo/RunParameters parsing updates, and SampleSheet v2 handling. Closes [#387](https://github.com/BU-ISCIII/iskylims/issues/387) [#392](https://github.com/BU-ISCIII/iskylims/pull/392)
- Improved Podman compatibility and adjusted SELinux bind mount handling in production compose setup [#393](https://github.com/BU-ISCIII/iskylims/pull/393)
- Fixed production container Django settings rendering when settings are provided through bind mounts.
- Fixed duplicate wetlab run configuration test logging.
- Fixed wetlab run configuration tests, including run discovery and completion checks.
- Fixed skipped wetlab run test states so they are displayed in the configuration test view.
- Fixed Samba cron path validation in wetlab configuration checks.
- Fixed Docker build cleanup so DNF cache cleanup only runs in the final DNF step.
- Improved error handling and logging when a SampleSheet cannot be copied to the documents directory.
- Fixed production container permissions for application logs, documents, static files, Django settings, temporary files, and Apache bind mounts.

#### Changed

- Updated installation script with variable modules for more flexibility [#269](https://github.com/BU-ISCIII/iskylims/pull/269)
- Updated installation script to remove commas in values of `rawtobunbarcode` table [#277](https://github.com/BU-ISCIII/iskylims/pull/277)
- Updated installation documentation and script, fixing small issues [#284](https://github.com/BU-ISCIII/iskylims/pull/284)
- Unify main and develop branches [#334](https://github.com/BU-ISCIII/iskylims/pull/334)
- Increased max upload memory size. (#328)
- Renamed method `get_delivery_date` to `get_delivered_date` for clarity.
- Improved query performance and excluded rejected/archived services from ongoing list. (#299)
- Updated sample metadata fields with standardized ontology mappings and schema alignment.[#358](https://github.com/BU-ISCIII/iskylims/pull/358)
- API create-sample-data Lab data mapping moved to core_config.LAB_REQUEST_ONTOLOGY_MAP [#377](https://github.com/BU-ISCIII/iskylims/pull/377)
- Renamed `docker-compose.yml` to `docker-compose.test.yml` and `docker-compose.prod.yml` for test clarity [#389](https://github.com/BU-ISCIII/iskylims/pull/389)
- Refactored Docker runtime handling when path is outside repository [#389](https://github.com/BU-ISCIII/iskylims/pull/389)
- Updated `docker_install.sh` to `container_install.sh` with multiple reliability improvements [#389](https://github.com/BU-ISCIII/iskylims/pull/389)
- Updated upgrade scripts documentation to include execution order information and docker upgrade clarifications [#389](https://github.com/BU-ISCIII/iskylims/pull/389)
- Improved wetlab API stats-info aggregation queries.
- Improved supercronic download resilience in the Docker build.
- Refreshed install fixtures during Docker upgrades.
- Updated production container documentation for bind mounts and generated configuration files.
- Normalized the Django settings bind path used by container installation.
- Updated installation configuration templates with clearer container configuration guidance.
- Simplified container runtime path handling by consistently using `INSTALL_PATH`.
- Updated production volume ownership handling for rootless Podman deployments.
- Updated the Spanish production guide and README with container permissions, bind mounts, configuration, and recovery procedures.

#### Removed

- Dummy fix in usage line [#271](https://github.com/BU-ISCIII/iskylims/pull/271)
- Removed migrations from `.gitignore` and app-level `.gitignore` to ensure version control of schema changes [#389](https://github.com/BU-ISCIII/iskylims/pull/389)
- Removed runtime ownership fixes from container installation.
- Removed the redundant `APP_INSTALL_PATH` container setting.

#### Requirements

| Package             | Last release Version  | New release Version  |
|:--------------------|:---------------------|:----------------------------|
| wheel               | 0.37.1               | 0.46.2                      |
| asn1crypto          | 1.5.0                | 1.5.1                       |
| bcrypt              | 4.0.1                | 4.2.0                       |
| biopython           | 1.79                 | 1.84                        |
| cryptography        | 38.0.3               | 44.0.3                      |
| Django              | 4.2                  | 4.2.28                      |
| django-crispy-forms | 2.0                  | 2.3                         |
| crispy-bootstrap5   |                      | 0.7                         |
| django-crontab      |                      | 0.7.1                       |
| django-js-asset     | 2.0.0                | 2.2.0                       |
| django-mptt         | 0.14.0               | 0.16.0                      |
| django-mptt-admin   | 2.4.1                | 2.6.2                       |
| django-cleanup      | 7.0.0                | 8.1.0                       |
| interop             |                      | >1.1.22                     |
| mod_wsgi            | 4.9.4                | 5.0.0                       |
| gunicorn            |                      | 22.0.0                      |
| mysqlclient         | 2.0.3                | 2.2.6                       |
| paramiko            | 3.1.0                | 3.4.1                       |
| jsonschema          | 4.17.3               | 4.23.0                      |
| pysmb               |                      | 1.2.9.1                     |
| django_extensions   | 3.2.1                | 3.2.3                       |
| djangorestframework | 3.14.0               | 3.15.2                      |
| drf-yasg            | 1.21.5               | 1.21.7                      |
| xlrd                |                      | 2.0.1                       |
| pandas              | 1.5.3                | 2.2.2                       |
| numpy               |                      | 1.26.4                      |
| openpyxl            | 3.1.1                | 3.1.5                       |
| setuptools          |                      | 78.1.1                      |
