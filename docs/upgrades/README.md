# Version-specific upgrade guides

Use these guides when an upgrade crosses a version boundary that requires special data migration or preparation. For upgrades without a listed special path, use the generic procedure in the main [README](../../README.md). Institutional rootless Podman deployments should also follow the current operational guidance in [LEAME](../../LEAME.md).

Always back up the database, application configuration, persistent documents, and installation files before upgrading.

## Upgrade path

Apply every relevant guide in order when skipping releases:

1. [2.3.0 to 3.0.0](2.3.0-to-3.0.0.md)
2. [3.0.0 to 3.1.0](3.0.0-to-3.1.0.md)
3. 3.1.0 to 3.1.1: no special data-migration scripts are required; use the generic upgrade procedure.

The supported historical starting point is written as `2.3.0`, rather than `2.x`, because that is the version explicitly documented by the migration scripts.

For a developer-oriented inventory of all migration scripts, see [UPGRADE_SCRIPTS.md](../../UPGRADE_SCRIPTS.md).
