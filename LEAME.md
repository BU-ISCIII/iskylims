# iSkyLIMS production installation and operations

## 1. Pre-installation checks

1. Record the reviewed revision for every service build context.
2. Verify host capacity, engine/Compose versions, DNS, TLS and network access.
3. Verify database and identity dependencies are reachable from containers.
4. Review the generated Compose model with the protected environment file.

## 2. Host preparation

Create only the documented database, document, log, configuration and proxy
bind sources. Apply the UID/GID, modes and SELinux labels declared by each
component's independent permission specification.

## 3. Configuration preparation

Prepare one mode-`0600` production settings file per application service.
Selected infrastructure add-ons contribute clearly marked sections to one
application's file; they do not introduce additional settings files.
Confirm that no file contains unresolved `CHANGE_ME` values and no `VITE_*`
setting contains a secret.

## 4. First installation

```bash
bash container_install.sh --action install --engine podman \
  --git_revision <reviewed-tag-or-commit> \
  --install_conf_map app,/protected/app_production_settings.txt
```

## 5. Verification

Run the generated smoke test against direct service endpoints and then verify
the public Apache/TLS endpoint, authentication redirects and one real read
workflow. Record `compose ps`, image IDs and deployed revisions.

## 6. Upgrade

1. Review the supported source/target versions and configuration differences.
2. Back up all assets listed below from one consistent recovery point.
3. Run `container_install.sh --action upgrade` with the new revision.
4. Stop on build, readiness, bootstrap, migration or smoke-test failure.

| Asset | Production location | Recovery requirement |
|---|---|---|
| `app` database | External production database | Database backup before migration |
| `app` documents | `app_documents` named volume | Volume backup |
| `app` static | `app_static` named volume | Replaceable through collectstatic |
| Samba test data | `samba_test_data` named volume | Disposable test/demo files |

## 7. Backup and restore

Test database dumps and filesystem/volume backups before relying on them.
Restore databases, documents and identity state from the same recovery point,
then deploy the compatible application revisions and rerun all smoke tests.

## 8. Rollback

Use an image-only rollback only when schema compatibility is confirmed.
Otherwise stop writes, restore the complete recovery point, deploy the recorded
previous revisions and verify direct and public endpoints before reopening.

## 9. Permission repair

```bash
bash container_install.sh --action fix-permissions --engine podman \
  --install_conf_map app,/protected/app_production_settings.txt
```

This action MUST NOT build images, migrate databases or delete data.

## 10. Incident troubleshooting

1. Capture service state and logs.
2. Test direct application health.
3. Test internal DNS and dependent databases/identity services.
4. Test Apache routing and external TLS last.
5. Preserve evidence and backups before any destructive recovery action.
