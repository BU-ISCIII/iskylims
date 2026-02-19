# Upgrade scripts

This file lists data migration scripts and the version range they apply to.
Run them with:

```bash
python manage.py runscript <script_name>
```

## 3.0.0 -> 3.1.0

- `wetlab/scripts/convert_rawtop_counter_to_int.py` (script name: `convert_rawtop_counter_to_int`)
- `wetlab/scripts/library_pool_to_many_relation.py` (script name: `library_pool_to_many_relation`)

## 2.3.0 -> 3.0.0

- `core/scripts/rename_app_name.py` (script name: `rename_app_name`)
- `core/scripts/migrate_sample_type.py` (script name: `migrate_sample_type`)
- `core/scripts/migrate_optional_values.py` (script name: `migrate_optional_values`)

## 2.3.0 -> 2.3.1

- `drylab/scripts/drylab_service_state_migration.py` (script name: `drylab_service_state_migration`)
