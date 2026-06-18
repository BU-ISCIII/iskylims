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
  - [Minimum requirements](#minimum-requirements)
  - [Docker deployment](#docker-deployment)
    - [Local test stack](#local-test-stack)
    - [Production container](#production-container)
      - [Persist logs/documents on the host](#persist-logsdocuments-on-the-host)
      - [Apache reverse proxy (container) + Gunicorn](#apache-reverse-proxy-container--gunicorn)
      - [Cron jobs inside the container](#cron-jobs-inside-the-container)
    - [Manage containers after installation](#manage-containers-after-installation)
    - [Upgrade docker deployment](#upgrade-docker-deployment)
    - [Upgrade docker deployment v3.0.0 to 3.1.0](#upgrade-docker-deployment-v300-to-310)
      - [Back up first](#back-up-first)
      - [Refresh code and settings](#refresh-code-and-settings)
  - [Bare-metal deployment (Ubuntu/CentOS)](#bare-metal-deployment-ubuntucentos)
    - [Install](#install)
      - [Clone the repository](#clone-the-repository)
      - [Prepare the database](#prepare-the-database)
      - [Configure install\_settings.txt](#configure-install_settingstxt)
      - [Run install.sh](#run-installsh)
    - [Upgrade (3.0.0 to 3.1.0)](#upgrade-300-to-310)
      - [Back up first](#back-up-first-1)
      - [Refresh code and settings](#refresh-code-and-settings-1)
      - [Run upgrade steps requiring root](#run-upgrade-steps-requiring-root)
      - [Run upgrade steps without root](#run-upgrade-steps-without-root)
  - [Common operations (Docker + bare-metal)](#common-operations-docker--bare-metal)
    - [Database creation, users and grants](#database-creation-users-and-grants)
    - [Backups](#backups)
    - [Restore / rollback](#restore--rollback)
    - [What to do if something fails](#what-to-do-if-something-fails)
  - [Final configuration steps](#final-configuration-steps)
    - [SAMBA configurarion](#samba-configurarion)
    - [Email verification](#email-verification)
  - [Developer notes](#developer-notes)
    - [Django migrations workflow](#django-migrations-workflow)
    - [Persistent host paths](#persistent-host-paths)
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

## Minimum requirements

Container deployment requirements:

- Docker Engine + Docker Compose v2, or Podman + `podman-compose`
- git >= 2.34 to clone/update the repository
- Host MySQL/MariaDB, Apache, Python, and `lsb_release` are not required for container deployment
- For local test containers: MySQL and Samba are started as containers by `container_install.sh --test`
- For production containers: access to an external MySQL/MariaDB server and Samba share configured in the selected install config
- Host directories and permissions for logs, documents, and static files, as described in [Persist logs/documents on the host](#persist-logsdocuments-on-the-host)

Bare-metal deployment requirements:

- **sudo privileges** for dependency installation
- MySQL >= 8.0 or MariaDB > 10.4
- Apache >= 2.4
- git >= 2.34
- Python >= 3.11
- Local email sender configured
- Access to the Samba share where run folders live
- `lsb_release` package:
  - RedHat/CentOS: `yum install redhat-lsb-core`
  - Ubuntu: `apt install lsb-core lsb-release`

## Docker deployment

### Local test stack

Bring up a full test stack (database, Samba, app) plus fixtures and demo data:

```bash
bash container_install.sh --test 2>&1 | tee test.log
```

Use `--engine podman` to run the same flow with Podman:

```bash
bash container_install.sh --test --engine podman 2>&1 | tee test.log
```

This uses `docker-compose.test.yml` by default.

Defaults can be customised:

- `--demo_data /path/to/iskylims_demo_data.tar.gz` to reuse a local demo archive (otherwise it is downloaded).
- `--skip_demo_data` or `--skip_test_data` to avoid loading extra data.
- `--install_type` (`full` by default) and `--git_revision` to control the build.
- `--script` to run one or more Django migration scripts through `install.sh` (repeat the flag as needed).

Example running a migration script during Docker install:

```bash
bash container_install.sh --test --script migrate_optional_values 2>&1 | tee test.log
```

When the script finishes, open `http://localhost:8001` and follow the prompt to create the Django superuser.

The image now includes the staged application tree under `${INSTALL_PATH}`. Test containers can therefore be recreated or restarted without rerunning the file installation step; only DB/bootstrap tasks are executed by `container_install.sh`.

### Production container

Deploy the iSkyLIMS container against external MySQL/Samba services:

1. Copy and edit the production settings template:

    ```bash
    cp conf/docker_production_settings.txt conf/my_prod_settings_iskylims.txt
    # edit conf/my_prod_settings_iskylims.txt with your DB/Samba details
    ```

2. Create the host directories used as bind mount sources and make them manageable by the account that will run `container_install.sh`.

   At minimum, this includes `/var/log/local/relecov-iskylims/apps`, `/var/log/local/relecov-iskylims/apache`, and any custom `APACHE_CONF_PATH` or `DJANGO_SETTINGS_PATH` parent directory set in `conf/my_prod_settings_iskylims.txt`.

    ```bash
    sudo mkdir -p /var/log/local/relecov-iskylims/apps /var/log/local/relecov-iskylims/apache
    sudo chown -R "$USER:$USER" /var/log/local/relecov-iskylims
    ```

   For rootless Podman, `container_install.sh` uses `podman unshare` to apply container UID/GID ownership where needed. For Docker, normal host permissions apply, so the script runner must be able to create and adjust the bind mount paths.

3. Build and run in production mode (uses `docker-compose.prod.yml` by default):

    ```bash
    bash container_install.sh --install_conf conf/my_prod_settings_iskylims.txt 2>&1 | tee ./iskylims_docker_install_$(date +%Y%m%d_%H%M%S).log
    ```

   Use `--compose_file` to override the compose file or `--install_type`/`--git_revision` to change the build.
   Add `--engine podman` to use Podman instead of Docker.
   Tip: capture logs for troubleshooting:

    ```bash
    bash container_install.sh --install_conf conf/my_prod_settings_iskylims.txt 2>&1 | tee ./iskylims_docker_install_$(date +%Y%m%d_%H%M%S).log
    ```

4. If this is a fresh install, create the Django superuser when prompted and complete the Samba configuration in the UI.

Production images now bake the staged iSkyLIMS application into the image itself. Host reboots or container recreation no longer require rerunning the app installation step; `container_install.sh` only performs runtime bootstrap tasks such as migrations, fixture refreshes, optional scripts, superuser creation on first install, and `collectstatic`.

Container build/runtime values are configured in the selected install config, not by exporting shell variables. Edit these fields in `conf/my_prod_settings_iskylims.txt` before running `container_install.sh`:

- `INSTALL_PATH`: runtime install root used by the app container, static/documents mounts, and install scripts. Default: `/opt/iskylims`.
- `APACHE_CONF_PATH`: host directory used for Apache bind-mounted config files. Leave empty to use `${INSTALL_PATH}/conf` as the host bind source; set this to a writable host path for rootless or hardened deployments.
- `SERVER_STATUS_SERVER_NAME`: Apache virtual host used for `/server-status`. Leave empty to use `DNS_URL`.
- `SERVER_STATUS_ALIASES`: aliases accepted by the server-status virtual host. Default: `127.0.0.1 localhost`.
- `SERVER_STATUS_ALLOW_FROM`: clients allowed to access `/server-status`. Default: `127.0.0.1 localhost`.
- `APACHE_FORWARDED_PROTO` / `APACHE_FORWARDED_PORT`: forwarded request scheme and port sent by Apache. Production defaults: `https` and `443`.
- `DJANGO_SETTINGS_PATH`: host path used for the bind-mounted Django `settings.py`. Leave empty to use `${INSTALL_PATH}/iskylims/settings.py` as the host bind source. If the value is a directory, ends with `/`, or does not end with `.py`, `container_install.sh` appends `settings.py`.
- `APP_UID` / `APP_GID`: runtime UID/GID for the `iskylims` user inside the container. Default: `1212:1212`.
- `APP_SHELL`: shell assigned to the runtime user during image build. Default: `/sbin/nologin`.
- `APP_PORT`: internal Gunicorn bind port for the `app` service. Default: `8001`.
- `DJANGO_DEBUG`: Django debug flag passed to the production app container. Default: `false`; keep it disabled in production.
- `DB_CONN_MAX_AGE`: Django persistent DB connection lifetime in seconds. Default: `60`.
- `WEB_CONCURRENCY`: Gunicorn worker count. Default: `2`.
- `GUNICORN_THREADS`: threads per Gunicorn worker. Default: `2`.
- `GUNICORN_TIMEOUT`: Gunicorn request timeout in seconds. Default: `300`.
- `GUNICORN_KEEPALIVE`: Gunicorn keep-alive in seconds. Default: `5`.

The standalone iSkyLIMS Apache container renders `APACHE_FORWARDED_PROTO` and `APACHE_FORWARDED_PORT` from this settings file. When iSkyLIMS runs inside the integrated RELECOV stack, the shared `relecov_apache` proxy uses the corresponding values from `my_prod_settings_relecov.txt`.

During production install/upgrade, `container_install.sh` writes `.env.prod.file` in the repository root. This file is ignored by git and is used by Compose for variable interpolation in `docker-compose.prod.yml`. It intentionally contains Compose/runtime metadata, not database or email passwords.

Host directory and ownership preparation is described in [Persist logs/documents on the host](#persist-logsdocuments-on-the-host).

#### Persist logs/documents on the host

The production compose file uses `INSTALL_PATH` from the selected install config as the app container runtime root. Apache config and Django settings bind sources use `APACHE_CONF_PATH` and `DJANGO_SETTINGS_PATH` when set.

Persistence layout:

- `/var/log/local/relecov-iskylims/apps` -> `${INSTALL_PATH}/logs` inside the `app` container
- `/var/log/local/relecov-iskylims/apache` -> `/var/log/httpd` inside the `apache` container
- `${APACHE_CONF_PATH:-${INSTALL_PATH}/conf}/iskylims_apache_reverse_proxy.conf` -> `/etc/httpd/conf.d/iskylims.conf` inside the `apache` container
- `${APACHE_CONF_PATH:-${INSTALL_PATH}/conf}/iskylims_apache_logs.conf` -> `/etc/httpd/conf.d/logformat.conf` inside the `apache` container
- `${APACHE_CONF_PATH:-${INSTALL_PATH}/conf}/iskylims_apache_server-status.conf` -> `/etc/httpd/conf.d/server-status.conf` inside the `apache` container
- `${DJANGO_SETTINGS_PATH:-${INSTALL_PATH}/iskylims/settings.py}` -> `${INSTALL_PATH}/iskylims/settings.py` inside the `app` container
- `iskylims_documents` named volume -> `${INSTALL_PATH}/documents`
- `iskylims_static` named volume -> `${INSTALL_PATH}/static`

If you override the compose file, ensure these mounts exist to keep logs and documents persistent.

Create host directories before the first deployment:

```bash
sudo mkdir -p /var/log/local/relecov-iskylims/apps
sudo mkdir -p /var/log/local/relecov-iskylims/apache
sudo mkdir -p <APACHE_CONF_PATH>
sudo mkdir -p <DJANGO_SETTINGS_PATH_PARENT>
sudo chown -R <APP_UID>:<APP_GID> /var/log/local/relecov-iskylims/apps <APACHE_CONF_PATH> <DJANGO_SETTINGS_PATH_PARENT>
```

For hardened/rootless Podman hosts, run the host preparation script as the same
user that starts the containers. The script pre-creates Apache log files, fixes
rootless Podman ownership for the UBI httpd user, and applies SELinux container
labels when SELinux is enabled:

```bash
bash hardening.sh
```

If an administrator runs it as root, set `PODMAN_USER` to the user that starts
the rootless containers:

```bash
PODMAN_USER=bioinfo bash hardening.sh
```

#### Apache reverse proxy (container) + Gunicorn

For production, the `app` container runs `gunicorn` (not `manage.py runserver`) and the `apache` service in `docker-compose.prod.yml` acts as the reverse proxy.

Static files:

- The app collects static files into `${INSTALL_PATH}/static`.
- `docker-compose.prod.yml` shares that directory with the `apache` service through the named volume `iskylims_static`.
- The reverse proxy config serves `/static` directly from `${INSTALL_PATH}/static`.

During `container_install.sh`, `conf/iskylims_apache_reverse_proxy.conf`, `conf/iskylims_apache_logs.conf`, and `conf/iskylims_apache_server-status.conf` are rendered and copied to `${APACHE_CONF_PATH}` on the host. If `APACHE_CONF_PATH` is empty, they are copied to `${INSTALL_PATH}/conf`. The reverse proxy `ServerName`, forwarded host, and access/error log file names are generated from `DNS_URL` in the selected install config. The `/server-status` virtual host uses `SERVER_STATUS_SERVER_NAME`, `SERVER_STATUS_ALIASES`, and `SERVER_STATUS_ALLOW_FROM`; by default it is restricted to localhost.

`container_install.sh` prepares a host-side Django `settings.py` bind source at `${DJANGO_SETTINGS_PATH}`, or at `${INSTALL_PATH}/iskylims/settings.py` when `DJANGO_SETTINGS_PATH` is empty. During the bootstrap step, `install.sh` updates that bind-mounted file from `conf/template_settings.txt` and the selected install config, preserving an existing `SECRET_KEY`. Runtime settings can then be edited and the container restarted without rebuilding the image.

If you need a different app container runtime root, set `INSTALL_PATH` in the install config file before running `container_install.sh`. If the runtime root is not writable on the host, set `APACHE_CONF_PATH` and `DJANGO_SETTINGS_PATH` to writable host paths in the install config file.

`container_install.sh` creates `${APACHE_CONF_PATH:-${INSTALL_PATH}/conf}`, the parent directory for `${DJANGO_SETTINGS_PATH:-${INSTALL_PATH}/iskylims/settings.py}`, `/var/log/local/relecov-iskylims/apps`, and `/var/log/local/relecov-iskylims/apache` before `compose up`, copies the three Apache config files there, prepares the bind-mounted Django settings file if it does not exist, passes runtime values into Compose, and then runs `install.sh --bootstrap ...` inside the `app` container. The container image already contains the staged Django project and virtualenv under `${INSTALL_PATH}`; the bootstrap step updates settings, applies migrations, optional scripts/fixtures, and refreshes `${INSTALL_PATH}/static`, while the Apache container keeps using the host log path `/var/log/local/relecov-iskylims/apache`.

SELinux note for pre-production and production:

- Ensure `/var/log/local/relecov-iskylims/apache` is writable by the container runtime and labeled for containers, for example `container_file_t`.
- If the host path is already labeled `container_file_t`, do not add `:Z` to the Apache log bind mount. `:Z` forces a relabel and may fail with `lsetxattr(... container_file_t ...): operation not permitted`.
- A quick check is:

```bash
ls -ldZ /var/log/local/relecov-iskylims/apache
```

- Expected example:

```text
system_u:object_r:container_file_t:s0
```

- If Apache fails on startup with `ModSecurity: Failed to open debug log file: /var/log/httpd/modsec_debug.log`, remove any stale host file and recreate/restart the container. In practice, deleting `/var/log/local/relecov-iskylims/apache/modsec_debug.log` has been enough when the existing inode had bad permissions/label state.

#### Cron jobs inside the container

Cron runs via `supercronic`, started by the container entrypoint script. The script reads Django `CRONJOBS` directly, writes `${INSTALL_PATH}/cron/iskylims`, and starts `supercronic` as the non-root app user. It does not call `manage.py crontab` or the system `crontab` command during container startup.

If you change `CRONJOBS`, rebuild or restart the container to regenerate the cron file.

### Manage containers after installation

After a production install, use the generated `.env.prod.file` whenever you run Compose directly. This keeps paths, UID/GID, ports, and Gunicorn settings aligned with the install config.

Docker Compose examples:

```bash
docker compose --env-file .env.prod.file -f docker-compose.prod.yml ps
docker compose --env-file .env.prod.file -f docker-compose.prod.yml logs --tail 200 app
docker compose --env-file .env.prod.file -f docker-compose.prod.yml restart app
docker compose --env-file .env.prod.file -f docker-compose.prod.yml up -d
```

Podman Compose examples:

```bash
podman compose --env-file .env.prod.file -f docker-compose.prod.yml ps
podman compose --env-file .env.prod.file -f docker-compose.prod.yml logs --tail 200 app
podman compose --env-file .env.prod.file -f docker-compose.prod.yml restart app
podman compose --env-file .env.prod.file -f docker-compose.prod.yml up -d
```

If you edit container runtime values in the install config, rerun `container_install.sh --install_conf <file>` so `.env.prod.file` and the running containers are regenerated consistently.

If containers were recreated manually, persistent volumes were restored, bind mount ownership changed, or `APP_UID` / `APP_GID` changed, repair permissions before running bootstrap tasks:

```bash
bash container_install.sh --engine podman --install_conf conf/my_prod_settings_iskylims.txt --action fix-permissions
```

This action does not rebuild images or run migrations. It refreshes `.env.prod.file` and fixes host bind mount permissions with `podman unshare` when needed. If `iskylims_app` is already running, it also fixes mounted app volumes from inside the container as root; otherwise, start the containers and rerun the same command to repair named volumes. For Docker, use `--engine docker`.

### Upgrade docker deployment

Keep the same `APP_UID`/`APP_GID` values in the selected install config before running an upgrade.

Re-deploy the application container against an existing production database:

```bash
bash container_install.sh --install_conf conf/my_prod_settings_iskylims.txt --action upgrade 2>&1 | tee ./iskylims_docker_install_$(date +%Y%m%d_%H%M%S).log
```

The upgrade path rebuilds/restarts the container and runs `install.sh --bootstrap upgrade --tables` inside the app container. The app files are already baked into the rebuilt image; the bootstrap phase applies migrations with `--fake-initial`, refreshes `conf/first_install_tables.json`, refreshes static files, and skips superuser/demo/test data loading.

### Upgrade docker deployment v3.0.0 to 3.1.0

#### Back up first

Run the backup steps in [Backups](#backups) first.

For 3.0.0 -> 3.1.0, export the LibraryPool mapping first, then run the upgrade with pre/post scripts:

```bash
mysql --user=<db_user> --password=<db_password> --host=<db_server_ip> --port=<db_port> iskylims \
  -e "SELECT id, run_process_id_id FROM wetlab_library_pool" \
  > /tmp/library_pool_run_process.tsv
```

#### Refresh code and settings

```bash
cd <your working directory>/iskylims
git pull
cp conf/docker_production_settings.txt my_prod_settings_iskylims.txt
sudo nano my_prod_settings_iskylims.txt
```

Ensure the file uses Linux-friendly encoding (UTF-8/ASCII) if you edit it on Windows.

Keep the same `APP_UID`/`APP_GID` values in the selected install config before running the 3.0.0 -> 3.1.0 upgrade.

Run upgrade command:

```bash
bash container_install.sh --engine podman --install_conf my_prod_settings_iskylims.txt --action upgrade \
  --script_before convert_rawtop_counter_to_int \
  --script_after library_pool_to_many_relation,/tmp/library_pool_run_process.tsv 2>&1 | tee ./iskylims_docker_install_$(date +%Y%m%d_%H%M%S).log
```

## Bare-metal deployment (Ubuntu/CentOS)

### Install

#### Clone the repository

```bash
cd <your working directory>
git clone https://github.com/BU-ISCIII/iskylims.git iskylims
cd iskylims
```

#### Prepare the database

Create the database and application user following [Database creation, users and grants](#database-creation-users-and-grants), then note DB host/port/user/password for `install_settings.txt`.

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

The staged/bootstrap split added for container images is internal. Bare-metal commands do not change: `--install` and `--upgrade` still run the complete dependency, application, and database workflow documented here.

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
  sudo bash install.sh --install full --git_revision main --tables 2>&1 | tee ./iskylims_install_$(date +%Y%m%d_%H%M%S).log
  ```

- If Apache is managed elsewhere, skip the automatic restart with `--skip_apache_restart`.

### Upgrade (3.0.0 to 3.1.0)

Follow these steps to move from version 3.0.0 to the 3.1.x series.

#### Back up first

Run the backup steps in [Backups](#backups) first.
- Additionally, back up the full installation folder (for example `/opt/iskylims`) for bare-metal rollback.
- If you use library pools, export them before upgrading:

  ```bash
mysql --user=<db_user> --password=<db_password> --host=<db_server_ip> --port=<db_port> iskylims \
  -e "SELECT id, run_process_id_id FROM wetlab_library_pool" \
  > /tmp/library_pool_run_process.tsv
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
bash install.sh --upgrade app --git_revision main \
  --script_before convert_rawtop_counter_to_int \
  --script_after library_pool_to_many_relation,/tmp/library_pool_run_process.tsv
```

Upgrades regenerate migrations and apply them with `--fake-initial` so existing tables remain intact, matching the Docker workflow.

## Common operations (Docker + bare-metal)

### Database creation, users and grants

Run as MySQL root:

```sql
CREATE DATABASE IF NOT EXISTS iskylims CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci;

CREATE USER IF NOT EXISTS 'iskylims'@'%' IDENTIFIED BY 'djangopass';
CREATE USER IF NOT EXISTS 'iskylims'@'localhost' IDENTIFIED BY 'djangopass';

GRANT ALL PRIVILEGES ON iskylims.* TO 'iskylims'@'%';
GRANT ALL PRIVILEGES ON iskylims.* TO 'iskylims'@'localhost';

FLUSH PRIVILEGES;
```

Verification:

```sql
SHOW GRANTS FOR 'iskylims'@'%';
```

### Backups

Database dump:

```bash
mysqldump -h <db_host> -P <db_port> -u iskylims -p iskylims > iskylims_$(date +%Y%m%d_%H%M%S).sql
```

Logs archive:

```bash
tar -czf iskylims_app_logs_$(date +%Y%m%d_%H%M%S).tgz -C /var/log/local/relecov-iskylims/apps .

tar -czf iskylims_apache_logs_$(date +%Y%m%d_%H%M%S).tgz -C /var/log/local/relecov-iskylims/apache .
```

Documents volume archive:

```bash
docker run --rm -v iskylims_documents:/from -v "$PWD":/to alpine \
  tar -czf /to/iskylims_documents_$(date +%Y%m%d_%H%M%S).tgz -C /from .
```

With Podman, use the same command replacing `docker` with `podman`.

Suggested order before upgrades:

1. DB dump
2. Documents volume archive
3. Logs archive

### Restore / rollback

Restore DB:

```bash
mysql -h <db_host> -P <db_port> -u iskylims -p iskylims < iskylims_YYYYMMDD_HHMMSS.sql
```

Restore documents volume:

```bash
docker run --rm -v iskylims_documents:/to -v "$PWD":/from alpine \
  sh -lc "cd /to && tar -xzf /from/iskylims_documents_YYYYMMDD_HHMMSS.tgz"
```

With Podman, use the same command replacing `docker` with `podman`.

Restore logs:

```bash
mkdir -p /var/log/local/relecov-iskylims/apps
tar -xzf iskylims_app_logs_YYYYMMDD_HHMMSS.tgz -C /var/log/local/relecov-iskylims/apps

mkdir -p /var/log/local/relecov-iskylims/apache
tar -xzf iskylims_apache_logs_YYYYMMDD_HHMMSS.tgz -C /var/log/local/relecov-iskylims/apache
```

Bare-metal full rollback example:

```bash
sudo rm -rf /opt/iskylims
sudo cp -r /home/dadmin/backup_prod/iSkyLIMS/ /opt/
sudo /scripts/hardening.sh
mysql -u iskylims -p -h <db_host> iskylims < /home/dadmin/backup_prod/bk_iSkyLIMS_YYYYMMDDHHMM.sql
```

### What to do if something fails

When install/upgrade fails, restore the previous state and retry with logs enabled.

Quick diagnostics:

```bash
# bare-metal
cd /opt/iskylims
python manage.py check

# docker
docker compose --env-file .env.prod.file -f docker-compose.prod.yml ps
docker compose --env-file .env.prod.file -f docker-compose.prod.yml logs --tail 200 app
# podman
podman compose --env-file .env.prod.file -f docker-compose.prod.yml ps
podman compose --env-file .env.prod.file -f docker-compose.prod.yml logs --tail 200 app
```

If you suspect a corrupted image/build cache in Docker:

```bash
docker compose --env-file .env.prod.file -f docker-compose.prod.yml build --no-cache app
docker compose --env-file .env.prod.file -f docker-compose.prod.yml up -d --force-recreate app
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

## Developer notes

### Django migrations workflow

Migrations are committed to the repo. Do not run `makemigrations` during install/upgrade.

Baseline + upgrade flow for new releases:

1. Generate baseline migrations from the last stable tag (example 3.0.0).
2. Commit the baseline migrations.
3. Generate new migrations on `develop` for schema changes and commit them.
4. Upgrades run `migrate --fake-initial` once to align existing tables, then `migrate` to apply the new migration files.

### Persistent host paths

See [Persist logs/documents on the host](#persist-logsdocuments-on-the-host) in the production deployment section.

### Configure Apache server

These steps apply to bare-metal Apache installations. Docker production deployments use the `apache` container described above and do not require copying configs into `/etc/apache2` or `/etc/httpd`.

Copy the apache configuration file according to your distribution inside the apache configuration directory and rename it to `iskylims.conf`.

Typical config locations:

- Ubuntu/Debian: `/etc/apache2/sites-available/iskylims.conf` (enable with `a2ensite`)
- CentOS/RHEL: `/etc/httpd/conf.d/iskylims.conf`

Suggested steps (host Apache as reverse proxy):

1. Copy the example config:

    ```bash
    sudo cp conf/iskylims_apache_reverse_proxy.conf /etc/apache2/sites-available/iskylims.conf
    # CentOS/RHEL:
    # sudo cp conf/iskylims_apache_reverse_proxy.conf /etc/httpd/conf.d/iskylims.conf
    ```

2. Edit the config:

    - Set `ServerName`
    - Ensure `ProxyPass` points to `http://localhost:8001/`
    - Ensure `Alias /static/ /opt/iskylims/static/`

3. Create the static folder on the host:

    ```bash
    sudo mkdir -p /opt/iskylims/static
    ```

4. Enable required modules (Ubuntu/Debian):

    ```bash
    sudo a2enmod proxy proxy_http headers
    sudo a2ensite iskylims.conf
    ```

5. Reload Apache:

    ```bash
    sudo systemctl reload apache2
    # CentOS/RHEL:
    # sudo systemctl reload httpd
    ```

### Verification of the installation

Open the navigator and type "localhost" or the "server local IP" and check that iSkyLIMs is running.

You can also check some of the functionality, while also checking samba and database connections using:

- Go to [configuration test](https://iskylims.isciii.es/wetlab/configurationTest/)
- Click submit
- Check all tabs so every connectin is successful.
- Run the 3 tests for each sequencing machine: MiSeq, NextSeq and NovaSeq.

## iSkyLIMS documentation

iSkyLIMS documentation is available at [https://iskylims.readthedocs.io/en/latest](https://iskylims.readthedocs.io/en/latest)
