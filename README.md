# iSkyLIMS

Laboratory information management system for wet-lab and dry-lab sequencing workflows.

## Architecture and services

| Service | Profile | Build context | Internal port |
|---|---|---|---:|
| `app` | `django` | `.` | settings: `APP_PORT` |

Profile behavior:

- Django services build with an ephemeral settings secret, render protected host settings, and run controlled migration/bootstrap steps.

Selected add-ons:

- Apache source configuration lives under `conf/apache/`; customize its virtual hosts and routes there. The installer renders final bind sources under `deployment/apache/`.
- The Samba add-on provides disposable NGS demo storage only in `--test` mode.

All selected services and add-ons are assembled into one
`docker-compose.prod.yml` and one `docker-compose.test.yml`.

## Requirements

- Docker with Compose v2, or Podman with a Compose provider.
- Git and access to every declared service build context.
- One protected production settings file for every application service.
- Production DNS, TLS termination, databases, storage and backup destinations
  required by the selected profiles.

## Configuration

Each application repository owns its configuration template. Copy it to a
protected, ignored file, set mode `0600`, and resolve every `CHANGE_ME` value.
Browser-visible `VITE_*` values MUST NOT contain secrets.
Apache and Keycloak settings appear as commented `APACHE_*` and `KEYCLOAK_*`
sections in the selected application's same settings file.

The installer generates one protected `.env.production.file` or
`.env.test.file`. Service prefixes prevent collisions between settings such as
`APP_PORT`, `DB_HOST`, `APP_UID` and `APP_GID`.

## Test installation

```bash
bash container_install.sh --test --action install --engine docker \
  --git_revision current
```

## Production installation

```bash
bash container_install.sh --action install --engine podman \
  --git_revision <reviewed-tag-or-commit> \
  --install_conf_map app,/protected/app_production_settings.txt
```

## Upgrade

Back up persistent state, review release/configuration changes, and run:

```bash
bash container_install.sh --action upgrade --engine podman \
  --git_revision <new-tag-or-commit> \
  --install_conf_map app,/protected/app_production_settings.txt
```

Django services run controlled bootstrap/migrations. React services rebuild an
immutable browser bundle and do not run Django bootstrap.

## Persistence and permissions

| Asset | Production location | Backup/rebuild policy |
|---|---|---|
| `app` database | External production database | Database backup before migration |
| `app` documents | `app_documents` named volume | Volume backup |
| `app` static | `app_static` named volume | Replaceable through collectstatic |
| Samba test data | `samba_test_data` named volume | Disposable test/demo files |

Run permission repair without building or migrating:

```bash
bash container_install.sh --action fix-permissions --engine podman \
  --install_conf_map app,/protected/app_production_settings.txt
```

Each application and add-on has independent host-bind and running-container
permission specifications in `container_install.sh`.

## Testing

```bash
bash scripts/smoke_test.sh --engine podman \
  --compose_file docker-compose.prod.yml \
  --env_file .env.production.file
```

The generated smoke dispatcher combines framework checks with direct HTTP
health checks. Add authenticated application workflows without removing the
generated baseline checks.

## Rollback and restore

Record the deployed image/revision before every upgrade. An image-only rollback
is safe only when the previous application version supports the current schema.
Otherwise restore databases and persistent files from the same recovery point.
See [LEAME.md](LEAME.md) for the ordered production procedure.

## Troubleshooting

- Inspect `compose ps` and service logs before restarting anything.
- Use `fix-permissions` for reviewed ownership/mode drift.
- Test direct service health before debugging Apache or external TLS.
- Verify Keycloak issuer, public hostname, clients and redirect URIs together.
- Never fake Django migrations or delete persistent volumes as a first response.
