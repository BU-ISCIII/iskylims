# Django schema migration workflow

This document defines how developers generate, consolidate, test, and release
Django migrations when the application is installed with
`container_install.sh`.

## Rules and terminology

- A **released migration** has been included in a release that may have been
  installed outside a disposable development environment. It is immutable.
- A **development migration** exists only on an unreleased branch. It may be
  replaced before release, but every database that applied it must then be
  rebuilt or restored to the last released migration state.
- A **release migration** is the reviewed migration delivered for one Django
  app in one release, for example `0002_release_2_0.py`.
- Migration numbers and names that have been released MUST NOT be reused.
- One migration per changed app per release is a target, not an absolute rule.
  Ordered data transformations, circular dependencies, staged deployments, or
  database limitations may require multiple files.
- `makemigrations` is a development command. Installation and upgrade MUST use
  committed migration files and MUST NOT generate migrations.

Django compares the current models with the state represented by migration
files. It does not derive migration operations by comparing models with the
live database.

## Installer actions and database state

`--test` selects disposable test infrastructure, but its named database volume
persists across container recreation. It does not mean that every invocation
starts with an empty database.

| Invocation | Required initial database state | Purpose |
|---|---|---|
| `--test --action install` | Empty test database | Test a new installation |
| `--test --action upgrade` | Recognized earlier migration state | Continue development or test an upgrade |
| Production `--action install` | Empty isolated or production database | New installation |
| Production `--action upgrade` | Recognized released migration state | Upgrade without discarding data |

Do not use `install` repeatedly against a populated database. Do not use
`upgrade` to adopt a legacy database with no migration history.

## Generate migrations without deploying the application

Generate migrations before running an installation. The normative portable
method uses a one-off image built from the current source tree; it does not
start the application service or connect to production.

The following example uses Docker. Podman accepts the same sequence with
`podman` in place of `docker`. Substitute the image name, install path, app
labels, and migration directories for the application.

Build a development image using the committed, non-secret test settings:

```bash
docker build \
  --build-arg INSTALL_CONF=conf/docker_test_settings.txt \
  --build-arg USE_INSTALL_CONF_SECRET=false \
  --build-arg RENDER_DJANGO_SETTINGS=true \
  --tag application-migrations:development .
```

Create a one-off container and generate meaningfully named migration files:

```bash
migration_container="$(docker create application-migrations:development \
  bash -lc 'cd "$INSTALL_PATH" && python manage.py makemigrations \
    --name release_2_0 app_one app_two')"
docker start --attach "$migration_container"
```

Copy only the generated migration directories back into the source checkout:

```bash
docker cp "$migration_container:/opt/application/app_one/migrations/." \
  app_one/migrations/
docker cp "$migration_container:/opt/application/app_two/migrations/." \
  app_two/migrations/
docker rm "$migration_container"
```

The `/opt/application` example MUST be replaced with the service's configured
`INSTALL_PATH`. Review ownership and permissions after copying files from a
rootful engine. Never copy generated settings, credentials, or the container's
complete staged tree into the source checkout.

Most projects can generate migrations without a running database. If custom
application imports or system checks require one, start only an isolated test
database dependency and point the one-off container at it. Migration generation
MUST NOT connect to production or require a deployed application service.

Review the new files, inspect their SQL where useful, and commit model changes
and their migrations together:

```bash
git diff -- app_one/models.py app_one/migrations \
  app_two/models.py app_two/migrations
git add app_one/models.py app_one/migrations \
  app_two/models.py app_two/migrations
git commit -m "Add schema migration for release 2.0"
```

Applications MAY document a Micromamba alternative if they provide all of the
following in the source checkout:

- A maintained environment definition containing the supported Django version
  and application dependencies.
- A source-tree `manage.py` and importable settings module.
- Development settings that cannot address production services or credentials.

Such an application may use, for example:

```bash
micromamba run -n application-development \
  python manage.py makemigrations --name release_2_0 app_one app_two
```

The standard does not assume that this shortcut exists. Generated Django
projects create `manage.py` and the deployment settings in the staged image, so
the one-off image method remains the portable workflow.

## Legacy database adoption

Use this procedure once when a supported database already contains application
tables but has no corresponding records in `django_migrations`.

1. Check out the last supported stable application version.
2. Generate and commit `0001_initial.py` from models matching that version.
3. Compare the complete schema represented by every initial migration with a
   restored production database. Check tables, columns, database types,
   nullability, defaults, primary and foreign keys, unique and check
   constraints, indexes, and many-to-many tables.
4. Back up the real database and retain the exact application revision used for
   verification.
5. Use a dedicated baseline release containing the verified initial migrations,
   but no schema changes for the next release.
6. Run `migrate --fake-initial` once through the documented adoption runbook.
7. Verify the recorded history and repeat the schema comparison.
8. Use normal container upgrades for every later release.

`--fake-initial` only performs limited object-existence checks. It does not
prove that types, constraints, indexes, defaults, or nullability match. It MUST
NOT replace the explicit comparison above and MUST NOT be used to hide a failed
or partially applied migration.

Legacy adoption is deliberately outside normal `container_install.sh install`
and `upgrade` operations. The application production runbook MUST document the
reviewed, application-specific command used inside the staged baseline image.

## First release

Before the first release, developers may create intermediate migrations while
models are changing. Because no released database depends on those names, the
files may be consolidated into a final `0001_initial.py` at release freeze.

1. Freeze model changes for the release candidate.
2. Preserve any deliberate `RunPython`, `RunSQL`, or custom operations.
3. Remove only unreleased migration files on the release branch.
4. Generate a new initial migration from the final models using the one-off
   development image.
5. Restore any required custom and data operations in the correct order.
6. Rebuild every disposable database that applied the deleted migrations.
7. Run the fresh-install container test described below.
8. Commit the models and `0001_initial.py` together.

Never regenerate an initial migration after it has been released.

## Consecutive releases

Assume release 1.0 contains:

```text
0001_initial.py
```

During development of 2.0, temporary files may accumulate:

```text
0001_initial.py
0002_add_status.py
0003_change_sample.py
0004_add_index.py
```

Developers with an existing test database apply new development migrations with:

```bash
bash container_install.sh --test --action upgrade --engine docker
```

At release freeze, consolidate only the unreleased files into a meaningful
release migration:

```text
0001_initial.py
0002_release_2_0.py
```

The released `0001_initial.py` remains byte-for-byte unchanged. After
temporarily removing only `0002` through `0004` from the release branch, use
the one-off image workflow above and run this command inside its container:

```bash
cd "$INSTALL_PATH"
python manage.py makemigrations --name release_2_0 app_one app_two
```

An approved project Micromamba environment may be used instead. Do not run the
command in the installer or an installed production application.

Consolidating final model state may omit operations that are not represented by
the models, including data copies, `RunPython`, `RunSQL`, temporary fields, and
carefully staged nullability changes. Review every intermediate migration and
carry required ordered operations into the release migration. Keep multiple
release migrations when that is safer.

Every database that applied a removed development migration is incompatible
with the consolidated names. Rebuild a disposable test database or restore a
production-like database to the last released snapshot before testing the final
migration. Never test the consolidated file by upgrading a database that still
records the deleted development migration names.

## Container release tests

Every release candidate MUST pass both a fresh installation and an upgrade from
the previous supported release. Use `container_install.sh` so tests exercise the
same build, bootstrap, migration, fixture, health, and smoke-test paths used by
operators.

### Fresh test installation

Start with confirmed-empty test database volumes and run:

```bash
bash container_install.sh --test --action install --engine docker
```

This proves that the complete committed migration graph creates the current
schema from an empty database. Recreating containers alone is insufficient
because named test volumes persist.

### Test upgrade

1. Install the previous released revision into an isolated test stack.
2. Load representative data or restore an anonymized production-like backup.
3. Preserve that test database volume.
4. Check out the release candidate.
5. Upgrade the existing stack:

```bash
bash container_install.sh --test --action upgrade --engine docker
```

Verify application data and constraints in addition to installer success.

### Production-path installation test

Use protected test credentials and an isolated empty database, then exercise
the production topology:

```bash
bash container_install.sh --action install --engine podman \
  --install_conf /protected/path/production-like-settings.txt
```

The database MUST be isolated from production even when it was restored from a
production backup.

### Production-path upgrade test

Restore the previous release's database into an isolated database, retain its
migration history and representative data, and run:

```bash
bash container_install.sh --action upgrade --engine podman \
  --install_conf /protected/path/production-like-settings.txt
```

Production deployment uses the same upgrade command after completing the
application runbook's backup, maintenance, and recovery-readiness steps.

## Release verification

Before tagging a release, verify all of the following:

- Model changes and migration files are committed together.
- `makemigrations --check --dry-run --noinput` reports no missing migrations.
- A fresh container installation succeeds from an empty database.
- A container upgrade succeeds from the previous supported release.
- `showmigrations --plan` contains no unapplied entries after each test.
- Data transformations have been exercised with representative data.
- Migration names recorded by tested databases exist in the committed graph or
  are valid predecessors declared by a reviewed squashed migration.
- Rollback or restore procedures have been tested for destructive operations.

## Long-term migration squashing

A long migration history is valid. When it becomes difficult to maintain,
Django's `squashmigrations` command may combine a reviewed range while declaring
which historical migrations it replaces. This is different from rewriting an
already released file.

1. Generate and review the squashed migration.
2. Keep the squashed migration and replaced files together for a compatibility
   release.
3. Test both a new installation and upgrades from every supported release.
4. Wait until supported installations have crossed the replaced history.
5. Remove the replaced files in a later release following Django's documented
   transition procedure.

Squashing is occasional maintenance, not part of every release. Migrations with
`RunPython`, `RunSQL`, custom operations, or cross-app dependencies require
additional review and may not squash safely.
