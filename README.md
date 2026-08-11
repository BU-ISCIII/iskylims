# iSkyLIMS

iSkyLIMS is a laboratory information management system for genomics facilities.
It tracks massive-sequencing work from library preparation and sequencing-run
registration through FASTQ generation, quality control, bioinformatics service
requests, analysis and delivery to researchers.

## What is iSkyLIMS?

iSkyLIMS connects wet-lab and dry-lab activities in one traceable workflow:

- **WetLab** records projects, samples, library preparation, pools and
  sequencing runs. It follows each run through registration, sample-sheet
  delivery, data processing, statistics generation and completion, and exposes
  reports by sample, project, run and investigation.
- **DryLab** manages bioinformatics service requests associated with sequencing
  data, including requested analyses, status, files, resolution and delivery.
- **Clinic** and shared core modules provide the supporting application data,
  permissions, configuration and APIs used by those workflows.

![iSkyLIMS sequencing workflow](img/iskylims_scheme.png)

## Infrastructure

The supported deployment separates application runtime from institutional
services and persistent data:

- The `app` service runs Django under Gunicorn as an unprivileged user.
- The `apache` service is the public container reverse proxy and serves the
  shared static files. Production TLS may terminate there or at the
  institution's upstream proxy, according to the reviewed deployment.
- Production uses an external MySQL or MariaDB database. The disposable test
  stack creates its own MySQL service.
- Production sequencing storage is an externally managed Samba share. The test
  stack provides a disposable Samba service and can load fixtures and demo NGS
  data through the installer.
- Documents and collected static files use persistent named volumes. Django
  settings, application logs and Apache logs use the standardized protected
  host paths documented below.
- Scheduled Django jobs run through Supercronic inside the application
  container from the project's `CRONJOBS` setting.

The same lifecycle supports Docker and Podman for local testing and production,
plus the reviewed Django bare-metal procedure. For issues, use the
[iSkyLIMS issue tracker](https://github.com/BU-ISCIII/iSkyLIMS/issues).

- [What is iSkyLIMS?](#what-is-iskylims)
- [Infrastructure](#infrastructure)
- [Get the code (required)](#get-the-code-required)
- [Choose your path](#choose-your-path)
- [Minimum requirements](#minimum-requirements)
- [Docker deployment](#docker-deployment)
  - [Local test stack](#local-test-stack)
  - [Production container](#production-container)
  - [Manage containers after installation](#manage-containers-after-installation)
  - [Upgrade docker deployment](#upgrade-docker-deployment)
- [Bare-metal deployment (Ubuntu/CentOS)](#bare-metal-deployment-ubuntucentos)
- [Common operations (Docker + bare-metal)](#common-operations-docker--bare-metal)
- [Final configuration steps](#final-configuration-steps)
  - [Configure Samba](#configure-samba)
  - [Verify email](#verify-email)
  - [Run the iSkyLIMS configuration tests](#run-the-iskylims-configuration-tests)
- [Developer notes](#developer-notes)
- [Application documentation](#application-documentation)

## Get the code (required)

```bash
git clone https://github.com/BU-ISCIII/iSkyLIMS.git relecov-iskylims
cd relecov-iskylims
```

For an orchestrated deployment, every external build context in the service
table must exist at the declared path relative to this checkout.

## Choose your path

| Capability | Supported | Owner or command |
|---|---:|---|
| Docker local test | Yes | `container_install.sh --test --engine docker` |
| Podman local test | Yes | `container_install.sh --test --engine podman` |
| Docker production | Yes | `container_install.sh --engine docker` |
| Podman production | Yes | `container_install.sh --engine podman` |
| Bare metal | Profile-specific | See [Bare-metal deployment](#bare-metal-deployment-ubuntucentos) |
| Upgrade | Yes | `--action upgrade` |
| Permission repair | Yes | `--action fix-permissions` |
| Backup and restore | Yes | Operator-owned; follow [LEAME.md](LEAME.md) |

Services:

| Service | Profile | Build context | Internal port |
|---|---|---|---:|
| `app` | `django` | `.` | settings: `APP_PORT` |

- Django services build with an ephemeral settings secret, render protected host settings, and run controlled migration/bootstrap steps.

Selected add-ons:

- Apache source configuration lives under `conf/apache/`; customize its virtual hosts and routes there. The installer renders final bind sources under `deployment/apache/`.
- The Samba add-on provides disposable NGS demo storage only in `--test` mode.

## Minimum requirements

- Git and access to every declared build context.
- Docker Engine with Compose v2, or Podman with a Compose provider.
- Enough disk and memory for image builds and persistent application data.
- A protected production settings file for every application service.
- Production DNS, TLS termination, database, storage, email, identity, backup,
  and monitoring services required by the selected profiles.

Copy each service's `conf/docker_production_settings.txt` to a protected,
ignored file, set mode `0600`, and replace every `CHANGE_ME` value. The exact
meaning and security classification of settings is in
[`conf/INSTALL_SETTINGS.md`](conf/INSTALL_SETTINGS.md).

## Docker deployment

Both engines use the same lifecycle and Compose files. Do not invoke Compose
directly for the first install or an upgrade: the installer also renders
configuration, prepares permissions, waits for readiness, and runs bootstrap.

### Local test stack

Docker:

```bash
bash container_install.sh --test --action install --engine docker \
  --git_revision current
```

Podman:

```bash
bash container_install.sh --test --action install --engine podman \
  --git_revision current
```

Test settings and test services are disposable. Verify either deployment with:

```bash
bash scripts/smoke_test.sh --test --engine docker
# or: bash scripts/smoke_test.sh --test --engine podman
```

Django test installation creates the disposable database declared by the test
Compose profile, waits for it, applies committed migrations, optionally loads
fixtures, runs selected data scripts, collects static files, and performs the
generated health checks.

Migration/data scripts are repeatable `django-extensions` runscript names. Use
`--script_before` for preparation before migrations and `--script` (an alias of
`--script_after`) for a transformation after migrations:

```bash
bash container_install.sh --test --action install --engine docker \
  --script_before prepare_test_data \
  --script migrate_optional_values
```

`--demo_data`, `--skip_demo_data`, and `--skip_test_data` are part of the
standard interface. A project that supplies fixtures or demo files must set
`application_supports_test_data=true` and implement `load_test_deployment_data`
in its wrapper; otherwise `--demo_data` is rejected explicitly.

For an automatic first administrator, set `CREATE_INITIAL_SUPERUSER=true` and
the `DJANGO_SUPERUSER_*` values in the selected test settings before install.
An existing account is never reset. Open the loopback URL using `APP_PORT` from
the rendered test environment, or `APACHE_PORT` when the Apache add-on is used.

The Samba add-on supplies disposable test storage only. Applications may load
fixtures and demo files into it through `load_test_deployment_data`; production
continues to use the externally managed storage configured by the application.

### Production container

Prepare the protected settings files and deploy a reviewed tag or commit.

Docker:

```bash
bash container_install.sh --action install --engine docker \
  --git_revision <reviewed-tag-or-commit> \
  --install_conf_map app,/protected/app_production_settings.txt
```

Podman:

```bash
bash container_install.sh --action install --engine podman \
  --git_revision <reviewed-tag-or-commit> \
  --install_conf_map app,/protected/app_production_settings.txt
```

The installer creates `.env.production.file`; use it for later direct Compose
operations. Production secrets stay in the protected settings files and are
not copied into image layers.

#### Persist logs/documents on the host

| Asset | Production location | Backup/rebuild policy |
|---|---|---|
| `app` database | External production database | Database backup before migration |
| `app` documents | `app_documents` named volume | Volume backup |
| `app` static | `app_static` named volume | Replaceable through collectstatic |
| `app` logs | `/var/log/local/relecov-iskylims/apps` host bind | Retain/rotate per institutional log policy |
| `app` rendered settings | `/srv/containers/bind/relecov-iskylims/settings/` host bind | Protected configuration backup |
| Apache logs | `/var/log/local/relecov-iskylims/apache` host bind | Retain/rotate per institutional log policy |
| Rendered Apache configuration | `deployment/apache/` in the deployment checkout | Rebuildable; preserve reviewed source configuration |
| Samba test data | `samba_test_data` named volume | Disposable test/demo files |

The standard fixes application binds below `/srv/containers/bind/relecov-iskylims`
and logs below `/var/log/local/relecov-iskylims`. The operator must still record the
backup owner, retention, actual engine volume names, and restore-test evidence
for every non-rebuildable asset. Never treat a container writable layer as
persistent storage.

#### Reverse proxy and application server

The selected profiles and add-ons define the internal application server and
proxy topology. Review public hostnames, TLS ownership, forwarded headers,
request limits, timeouts, health paths, and static/media routing together.

#### Scheduled jobs

The application developer must list every scheduler/worker, whether a failed
job blocks a workflow, and how operators inspect and retry it. Do not add an
untracked host cron job when the application profile owns scheduling.

### Manage containers after installation

Use the engine that performed the installation:

```bash
docker compose --env-file .env.production.file -f docker-compose.prod.yml ps
docker compose --env-file .env.production.file -f docker-compose.prod.yml logs --tail 200
docker compose --env-file .env.production.file -f docker-compose.prod.yml restart
```

```bash
podman compose --env-file .env.production.file -f docker-compose.prod.yml ps
podman compose --env-file .env.production.file -f docker-compose.prod.yml logs --tail 200
podman compose --env-file .env.production.file -f docker-compose.prod.yml restart
```

### Upgrade docker deployment

After taking a consistent backup and reading the version-specific upgrade
notes:

```bash
bash container_install.sh --action upgrade --engine podman \
  --git_revision <new-reviewed-tag-or-commit> \
  --install_conf_map app,/protected/app_production_settings.txt
```

Replace `podman` with `docker` for a Docker-managed deployment. Stop on build,
readiness, bootstrap, migration, or smoke-test failure. See [LEAME.md](LEAME.md)
for the ordered production checklist and rollback decision.

## Bare-metal deployment (Ubuntu/CentOS)

### Install

#### Clone the repository

Use [Get the code (required)](#get-the-code-required) and check out the reviewed
revision.

#### Prepare the database

Provision the application database and least-privilege account outside the
installer. Confirm that the host can reach it before bootstrap.

#### Configure install_settings.txt

Start from `conf/docker_production_settings.txt`, but review all paths and
container-oriented defaults for the target host. Keep the resulting file
ignored and mode `0600`.

#### Run install.sh

The Django profile includes `install.sh` for application staging and bootstrap,
but system package, database, web-server, service-manager, TLS, and backup
provisioning remain host-specific. Bare-metal installation is supported only
after the application developer documents and tests those integrations.

```bash
# Stage application files and dependencies.
bash install.sh --stage install --git_revision current \
  --conf conf/docker_production_settings.txt

# Bootstrap the prepared runtime (settings, migrations and static files).
bash install.sh --bootstrap install \
  --conf conf/docker_production_settings.txt
```

For upgrades, take a backup and replace both `install` actions with `upgrade`.
Do not use container-oriented paths or defaults on a bare-metal host without an
application-specific review.

For a host-managed Apache 2.4 deployment, adapt the reviewed virtual host from
`conf/apache/` to the distribution path. The generated add-on files target the
container image, so do not copy them blindly without checking module names,
paths, runtime user, TLS ownership, and log locations.

Ubuntu/Debian baseline:

```bash
sudo cp <reviewed-apache-vhost.conf> /etc/apache2/sites-available/relecov-iskylims.conf
sudo a2enmod proxy proxy_http headers
sudo a2ensite relecov-iskylims.conf
sudo apache2ctl configtest
sudo systemctl reload apache2
```

CentOS/RHEL baseline:

```bash
sudo cp <reviewed-apache-vhost.conf> /etc/httpd/conf.d/relecov-iskylims.conf
sudo httpd -t
sudo systemctl reload httpd
```

The reviewed virtual host must define the public `ServerName`, proxy to the
Django `APP_PORT`, serve the correct static/media paths, preserve forwarded
scheme/host headers, and use institutionally managed TLS and logs.

### Upgrade bare-metal deployment

Follow the same staged lifecycle with `upgrade` only after a consistent backup
and review of the version-specific guide.

## Common operations (Docker + bare-metal)

### Database creation, users and grants

Production databases are externally managed unless the application documents a
different supported topology. Create a dedicated schema and least-privilege
account, verify connectivity from the application container, and keep DBA
commands and credentials outside this repository.

Connect as an authorized database administrator without putting the password
on the command line:

```bash
DB_HOST='CHANGE_ME'
DB_PORT='3306'
DB_ADMIN='CHANGE_ME'
DB_NAME='CHANGE_ME'
DB_USER='CHANGE_ME'
mysql --host="$DB_HOST" --port="$DB_PORT" --user="$DB_ADMIN" --password
```

Create the application database and account. Replace every angle-bracket value;
restrict the account host further than `%` when the network topology permits.

```sql
CREATE DATABASE `<db-name>`
  CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci;
CREATE USER '<db-user>'@'%' IDENTIFIED BY '<strong-generated-password>';
GRANT ALL PRIVILEGES ON `<db-name>`.* TO '<db-user>'@'%';
FLUSH PRIVILEGES;
```

Verify the same endpoint and least-privilege credentials configured for the
application:

```bash
mysql --host="$DB_HOST" --port="$DB_PORT" --user="$DB_USER" --password \
  --database="$DB_NAME" --execute='SELECT 1;'
```

### Backups

Back up every non-rebuildable row in the persistence table from one consistent
recovery point before installation or upgrade. Record the revision, image IDs,
settings files, and backup identifiers.

```bash
BACKUP_DIR="/srv/containers/backup/relecov-iskylims/$(date +%Y%m%d_%H%M%S)"
SETTINGS_FILE='/protected/app_production_settings.txt'
DOCUMENTS_VOLUME='CHANGE_ME'
DB_HOST='CHANGE_ME'
DB_PORT='3306'
DB_NAME='CHANGE_ME'
DB_USER='CHANGE_ME'
mkdir -p "$BACKUP_DIR"
git rev-parse HEAD > "$BACKUP_DIR/git-revision.txt"
cp .env.production.file "$SETTINGS_FILE" "$BACKUP_DIR/"
chmod -R go-rwx "$BACKUP_DIR"

mysqldump --single-transaction --routines --triggers \
  --host="$DB_HOST" --port="$DB_PORT" --user="$DB_USER" --password \
  "$DB_NAME" > "$BACKUP_DIR/database.sql"

podman volume ls | grep 'relecov-iskylims'
podman volume export "$DOCUMENTS_VOLUME" > "$BACKUP_DIR/documents.tar"
tar -C /srv/containers/bind -czf "$BACKUP_DIR/bind-mounts.tar.gz" relecov-iskylims
sha256sum "$BACKUP_DIR"/* > "$BACKUP_DIR/SHA256SUMS"
```

For Docker, archive a named volume through a temporary container after ensuring
the application is not writing to it:

```bash
docker run --rm \
  --volume "$DOCUMENTS_VOLUME":/data:ro \
  --volume "$BACKUP_DIR":/backup \
  alpine tar -C /data -cf /backup/documents.tar .
```

The full ordered backup checklist, including logs and image metadata, is in
[LEAME.md](LEAME.md).

### Restore / rollback

An image-only rollback is safe only when the previous application version
supports the current schema and persistent-file format. Otherwise stop writes,
restore the database and files from the same recovery point, deploy the recorded
compatible revision, and rerun all smoke tests.

Compatible application-only rollback:

```bash
bash container_install.sh --action upgrade --engine podman \
  --git_revision <previous-reviewed-revision> \
  --install_conf_map app,/protected/app_production_settings.txt
```

Full restore when schema or persistent-file formats are incompatible:

```bash
BACKUP_DIR='/srv/containers/backup/relecov-iskylims/CHANGE_ME'
DOCUMENTS_VOLUME='CHANGE_ME'
DB_HOST='CHANGE_ME'
DB_PORT='3306'
DB_NAME='CHANGE_ME'
DB_USER='CHANGE_ME'
podman compose --env-file .env.production.file -f docker-compose.prod.yml down
mysql --host="$DB_HOST" --port="$DB_PORT" --user="$DB_USER" --password \
  "$DB_NAME" < "$BACKUP_DIR/database.sql"
podman volume import "$DOCUMENTS_VOLUME" "$BACKUP_DIR/documents.tar"
tar -C /srv/containers/bind -xzf "$BACKUP_DIR/bind-mounts.tar.gz"
bash container_install.sh --action fix-permissions --engine podman \
  --install_conf_map app,/protected/app_production_settings.txt
```

Then deploy the revision recorded in `git-revision.txt`, start the deployment,
and run the smoke test before reopening service. For Docker volume restoration,
reverse the temporary-container archive command by mounting the empty target
volume at `/data` and extracting `/backup/documents.tar` there.

### What to do if something fails

1. Preserve installer output, `compose ps`, image IDs, and service logs.
2. Test the direct application health endpoint and dependencies.
3. Test proxy routing, public DNS, and TLS after direct health succeeds.
4. Run permission repair for reviewed ownership or SELinux drift:

   ```bash
   bash container_install.sh --action fix-permissions --engine podman \
     --install_conf_map app,/protected/app_production_settings.txt
   ```

5. Do not fake migrations, delete volumes, or rebuild from an unrecorded
   revision as a first response.

### Service-specific operational commands

#### Django service `app`

```bash
# Logs and an interactive shell (replace podman with docker when applicable).
podman compose --env-file .env.production.file -f docker-compose.prod.yml \
  logs --tail 200 app
podman compose --env-file .env.production.file -f docker-compose.prod.yml \
  exec app bash

# Rebuild static assets without running migrations.
podman compose --env-file .env.production.file -f docker-compose.prod.yml \
  exec app bash -lc \
  'cd "$INSTALL_PATH" && source virtualenv/bin/activate && python manage.py collectstatic --noinput'

# Inspect Django and migration state before deciding whether to recover.
podman compose --env-file .env.production.file -f docker-compose.prod.yml \
  exec app bash -lc \
  'cd "$INSTALL_PATH" && source virtualenv/bin/activate && python manage.py check --deploy && python manage.py showmigrations --plan'
```

For bootstrap recovery, fix the cause and rerun `container_install.sh` with the
same revision, protected configuration, and `--action install` or `upgrade`.
This safely recreates the temporary runtime configuration and repeats the
controlled migration/fixture/static lifecycle. Direct `manage.py migrate` is a
diagnostic last resort and must use the same backup and release procedure.

#### Apache service

```bash
podman compose --env-file .env.production.file -f docker-compose.prod.yml \
  logs --tail 200 apache
podman compose --env-file .env.production.file -f docker-compose.prod.yml \
  exec apache httpd -t

APACHE_PORT='CHANGE_ME'
SERVER_STATUS_SERVER_NAME='localhost'
curl --fail --show-error \
  --header "Host: $SERVER_STATUS_SERVER_NAME" \
  "http://127.0.0.1:$APACHE_PORT/server-status?auto"
```

Keep `SERVER_STATUS_ALLOW_FROM` limited to trusted diagnostic hosts. If SELinux
is enabled, inspect the persistent log bind and confirm a container-compatible
label before restarting:

```bash
ls -ldZ /var/log/local/relecov-iskylims/apache
```

An Apache failure containing `ModSecurity: Failed to open debug log file` often
means the existing `modsec_debug.log` inode has stale ownership or labeling.
Preserve it for diagnosis, run `fix-permissions`, and restart Apache. If it must
be replaced, move it to a timestamped backup instead of deleting evidence:

```bash
sudo mv /var/log/local/relecov-iskylims/apache/modsec_debug.log \
  /var/log/local/relecov-iskylims/apache/modsec_debug.log.blocked
bash container_install.sh --action fix-permissions --engine podman \
  --install_conf_map app,/protected/app_production_settings.txt
podman compose --env-file .env.production.file -f docker-compose.prod.yml restart apache
```

## Final configuration steps

Sign in with the administrator account after the first installation. The
installer creates that account only when `CREATE_INITIAL_SUPERUSER=true` and
the protected `DJANGO_SUPERUSER_*` settings are provided.

### Configure Samba

iSkyLIMS stores the connection to the institutional sequencing-data share in
its application configuration. The Samba container included in the test stack
is disposable; production must use the reviewed external share.

1. Open **Massive sequencing** from the iSkyLIMS home page.

   ![Open the WetLab module](img/got_to_wetlab.png)

2. Go to **Configuration → Samba configuration**.
3. Enter the server/share, credentials and paths for the sequencing run
   folders. Use the disposable `SAMBA_USER` and `SAMBA_PASSWORD` only in the
   local test deployment.
4. Save the configuration and verify the connection before enabling scheduled
   run discovery.

   ![iSkyLIMS Samba configuration form](img/samba_form.png)

### Verify email

1. Open **Massive sequencing → Configuration → Email configuration**.
2. Confirm the sender and recipients expected by the laboratory workflow.
3. Send a test message and verify delivery. Do not approve production while
   application email fails, even when the SMTP connection itself succeeds.

### Run the iSkyLIMS configuration tests

After every clean installation, upgrade, restore or Samba/email change:

1. Run the generated deployment smoke test:

   ```bash
   bash scripts/smoke_test.sh --engine podman
   ```

   Replace `podman` with `docker` for a Docker-managed deployment.

2. Sign in as an administrator and open
   `<public-iskylims-url>/wetlab/configurationTest/`. For the institutional
   deployment, the current endpoint is
   [iSkyLIMS configuration test](https://iskylims.isciii.es/wetlab/configurationTest/).
3. Submit the configuration test and inspect every result tab.
4. Confirm the database and Samba checks succeed.
5. Run the available checks for each configured sequencing platform, including
   MiSeq, NextSeq and NovaSeq where enabled.
6. Confirm a representative WetLab read workflow and record the result with the
   deployed revision.

## Developer notes

### Shared container installer library

`container_install.sh` sources the vendored files under
`deployment/lib/container/`. Do not edit those copies. Check or update them
from the standards repository with `scaffold.py check-lib` or `sync-lib`.

### Schema migration workflow

Django migrations MUST be generated, reviewed, tested, and committed with the
release. Installation and production upgrade run `migrate --noinput`; they
MUST NOT run `makemigrations` or silently manufacture schema history.

For a legacy application entering the standard:

1. Generate and commit baseline migrations from the last supported stable tag.
2. Generate and commit new migrations for later model changes.
3. Identify only legacy application labels whose existing tables match the
   committed initial migrations exactly.
4. Document and run `migrate <app-label> --fake-initial` once through a reviewed
   application migration callback or version-specific upgrade step. The common
   bootstrap does not apply `FAKEINITIAL_MODULES` automatically.
5. Put ordered data transformations in version-specific upgrade guides and run
   them through `--script_before`, `--script_after`, or `--script`.
6. Verify `showmigrations --plan` has no unapplied entries after bootstrap.

Never use `--fake` to conceal a failed or partially applied migration. Normal
new installations and subsequent upgrades use the committed migration graph
without `--fake-initial`.

### Persistent host paths

Keep source checkouts, protected configuration, bind mounts, engine-managed
volumes, logs, and backups separate. For rootless Podman, run the installer as
the same unprivileged account every time and use `fix-permissions` instead of
manually changing engine storage.

### Verification of the installation

```bash
bash scripts/smoke_test.sh --engine podman
```

Complete the application-specific Samba, email and sequencing checks in
[Final configuration steps](#final-configuration-steps) after the generated
baseline succeeds.

## Application documentation

- [User and administrator documentation](https://iskylims.readthedocs.io/en/latest/)
- [Version-specific upgrade guides](docs/upgrades/README.md)
- [Issue tracker](https://github.com/BU-ISCIII/iSkyLIMS/issues)
