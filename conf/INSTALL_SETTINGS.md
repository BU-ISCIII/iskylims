# Installation settings for iSkyLIMS

`docker_test_settings.txt` contains disposable local defaults.
`docker_production_settings.txt` is a template and MUST NOT contain real
production secrets. Operators copy it to an ignored, permission-restricted file.

Values already rendered from `project.json` (application name, module, Python
version and default paths) provide a runnable baseline. Operators customize
environment-dependent values. Application developers add domain settings to
both settings templates, this matrix, and `template_settings.py`; adding an
undocumented environment variable alone does not configure Django.

## Application and filesystem

| Variable | Required | Secret | Meaning |
|---|---:|---:|---|
| `REPO_PATH` | yes | no | Staged application source inside the image/container |
| `INSTALL_PATH` | yes | no | Application runtime root inside the container |
| `PROJECT_MODULE` | generated | no | Django package declared by `project.json`; do not customize independently |
| `PYTHON_BIN_PATH` | yes | no | Python used to create the virtual environment |
| `REQUIRED_MODULES` | application | no | Import checks required before bootstrap |
| `MIGRATION_MODULES` | application | no | Modules whose committed migrations are applied |
| `APP_UID`, `APP_GID` | yes | no | Runtime identity and rootless volume ownership |
| `APP_SHELL` | yes | no | Runtime account shell; normally `/sbin/nologin` in production |
| `APP_PORT` | yes | no | Internal Gunicorn and host-loopback port |
| `HOST_LOG_PATH` | production | no | Persistent application logs on the host |
| `DJANGO_SETTINGS_PATH` | production | no | Protected rendered `settings.py` bind source on the host |

## Database

| Variable | Required | Secret | Meaning |
|---|---:|---:|---|
| `DB_HOST`, `DB_PORT` | yes | no | MySQL endpoint reachable by the application |
| `DB_NAME`, `DB_USER` | yes | no | Application schema and least-privilege user |
| `DB_PASSWORD` | yes | yes | Application database password |
| `DB_ROOT_PASSWORD` | test only | yes | Root password for disposable Compose MySQL |

Production uses an external database. The production Compose file deliberately
contains no database service and publishes no database port. When that database
runs on the container host, use `host.docker.internal` with either Docker or
Podman; the production service maps it to the host gateway.

## Django and HTTP

| Variable | Required | Secret | Meaning |
|---|---:|---:|---|
| `DJANGO_DEBUG` | yes | no | MUST be `false` in production |
| `DJANGO_SECRET_KEY` | yes | yes | Stable Django signing key; preserve on upgrade |
| `DJANGO_ALLOWED_HOSTS` | yes | no | Exact production hostnames |
| `DJANGO_CSRF_TRUSTED_ORIGINS` | production | no | Exact HTTPS origins |
| `DB_CONN_MAX_AGE` | yes | no | Django persistent database connection lifetime |
| `WEB_CONCURRENCY` | production | no | Gunicorn workers |
| `GUNICORN_THREADS` | production | no | Threads per worker |
| `GUNICORN_TIMEOUT` | production | no | Request timeout in seconds |
| `GUNICORN_KEEPALIVE` | production | no | Keep-alive time in seconds |
| `APP_START_WAIT_TIMEOUT_SECONDS` | yes | no | Maximum wait for staged application files |
| `LOG_TYPE`, `LOG_PATH` | application | no | Application logging backend and optional location |

## Email

`EMAIL_HOST`, `EMAIL_PORT`, `EMAIL_HOST_USER`, `EMAIL_HOST_PASSWORD`, and
`EMAIL_USE_TLS` configure SMTP. The password is secret; the remaining values
are operational unless the username is sensitive locally. Email settings are
required when the application sends operational or account messages. Document
whether failed email blocks the user workflow and add an SMTP test to
production acceptance.

## Initial administrator

`CREATE_INITIAL_SUPERUSER=true` creates the first Django administrator only
during `--bootstrap install`. Set `DJANGO_SUPERUSER_USERNAME`, optional
`DJANGO_SUPERUSER_EMAIL`, and secret `DJANGO_SUPERUSER_PASSWORD` in the
protected settings file. Bootstrap retries never reset an existing account.

## Project-specific settings

Add every application setting here before declaring installation complete.
Classify secrets and identify cross-service values that must match an identity
provider, proxy, worker, or frontend.

Developer review checklist:

- map every project-specific token in `template_settings.py`;
- provide safe disposable test values and `CHANGE_ME` production placeholders;
- state validation rules, owner, restart/rebuild impact, and secret status;
- add acceptance checks for email, identity, storage, workers, and scheduled
  jobs used by real workflows.

## Selected infrastructure add-ons

Add-ons use independent settings below `conf/<addon>/`. For production, copy
the required add-on templates to protected files and pass them through the same
repeatable `--install_conf_map <component>,<path>` option used by application
services.

### Apache

`APACHE_LOG_PATH`, `APACHE_BIND_HOST`, `APACHE_PORT`,
`APACHE_FORWARDED_PROTO`, `APACHE_FORWARDED_PORT`, and
`APACHE_LIMIT_REQUEST_BODY` configure the rendered proxy. They are
operational values, not Django or React application settings.

`APACHE_SERVER_NAME` is the host name handled by the baseline VirtualHost.
`APACHE_UPSTREAM_SERVICE` defaults to `ADDONS.apache.CONFIG_SERVICE`, while
`APACHE_UPSTREAM_PORT` defaults to that service's `APP_PORT`.
`APACHE_PROXY_TIMEOUT` defaults to its `GUNICORN_TIMEOUT` (or 120 seconds), and
`APACHE_LOG_STEM` defaults to a filename-safe form of `APACHE_SERVER_NAME`.
Leave those four derived values empty unless the proxy route needs an override.

`SERVER_STATUS_SERVER_NAME`, `SERVER_STATUS_ALIASES`, and
`SERVER_STATUS_ALLOW_FROM` configure the restricted Apache status endpoint.
Keep its allow-list limited to trusted diagnostic hosts.

Edit the source files under `conf/apache/` to define the application's virtual
hosts, routes, and aliases. During installation they are rendered with the
protected deployment environment into `deployment/apache/`; only those final
files are bind-mounted. `APACHE_LOG_PATH` is the writable persistent host log
source.


### Samba test data

`SAMBA_USER` and `SAMBA_PASSWORD` configure the disposable Samba service used
only by the test Compose profile. The add-on creates no production service.
