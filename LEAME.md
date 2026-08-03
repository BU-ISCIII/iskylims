# Actualizacion de iSkyLIMS con Podman rootless

Esta guia es para la actualizacion del despliegue institucional de iSkyLIMS usando Podman rootless. Se asume que iSkyLIMS ya existe en la institucion, que la base de datos de produccion ya esta creada y que el objetivo es actualizar o recrear los contenedores de aplicacion sin crear una base de datos desde cero.

`container_install.sh` se encarga de construir la imagen, arrancar los contenedores, preparar permisos, actualizar configuracion, aplicar migraciones, refrescar estaticos y ejecutar los pasos necesarios de actualizacion.

## Indice

- [Actualizacion de iSkyLIMS con Podman rootless](#actualizacion-de-iskylims-con-podman-rootless)
  - [Indice](#indice)
  - [Requisitos](#requisitos)
  - [Estructura de directorios en los servidores](#estructura-de-directorios-en-los-servidores)
  - [Preparar directorios del host](#preparar-directorios-del-host)
  - [Actualizar codigo](#actualizar-codigo)
  - [Configurar `my_prod_settings_iskylims.txt`](#configurar-my_prod_settings_iskylimstxt)
  - [Backup antes de actualizar](#backup-antes-de-actualizar)
  - [Ejecutar la actualizacion](#ejecutar-la-actualizacion)
  - [Caso especial: actualizacion desde 3.0.0 a 3.1.0](#caso-especial-actualizacion-desde-300-a-310)
  - [Comprobaciones posteriores](#comprobaciones-posteriores)
  - [Rollback](#rollback)
  - [Reparar permisos](#reparar-permisos)
  - [Operaciones utiles](#operaciones-utiles)
  - [Notas de permisos](#notas-de-permisos)

## Requisitos

El despliegue usa:

- Podman rootless ejecutado por un usuario normal del sistema.
- `podman-compose` o `podman compose`.
- `container_install.sh` desde el repositorio de iSkyLIMS.
- `docker-compose.prod.yml`, lanzado con Podman.
- Una base de datos MySQL/MariaDB externa ya existente.
- Una configuracion Samba/storage ya existente o configurada desde la interfaz.
- Bind mounts del host para logs, configuracion Apache y `settings.py`.
- Volumenes Podman para `documents` y `static`.

Comprueba que Podman funciona sin root:

```bash
podman info
podman ps
```

No ejecutes `container_install.sh` con `sudo`. El usuario que ejecuta Podman debe ser el mismo usuario que ejecuta `container_install.sh`.

Los ejemplos usan `podman compose`. Si tu servidor solo tiene `podman-compose`, sustituye `podman compose` por `podman-compose`.

## Estructura de directorios en los servidores

Los servidores de desarrollo, preproduccion y produccion siguen la misma convencion de directorios. Cada tipo de informacion tiene una ubicacion concreta para separar el codigo y la configuracion del despliegue, los datos gestionados por Podman, los bind mounts y los logs. La estructura siguiente es orientativa: solo muestra una aplicacion como ejemplo y no pretende enumerar todos los directorios o ficheros existentes.

```text
/opt/containers_apps/
└── iskylims/
    ├── backup/                 # Opcional: backups propios del despliegue
    └── iskylims/               # Clon Git, fuentes y configuracion de instalacion

/srv/containers/
├── backup/                     # Backups centralizados, si no estan junto al despliegue
├── bind/
│   └── iskylims/               # Bind mounts organizados por aplicacion
│       ├── iskylims_apache_conf/
│       ├── iskylims_app_setting/
│       └── iskylims_django_setting/
├── shared/                     # Datos compartidos entre aplicaciones, cuando proceda
└── storage/
    └── <usuario-podman>/       # Almacenamiento interno rootless de Podman
        ├── overlay/
        ├── overlay-containers/
        ├── overlay-images/
        └── volumes/            # Volumenes persistentes gestionados por Podman

/var/log/local/
└── iskylims/
    ├── apache/                 # Logs del servidor web y de ModSecurity
    └── apps/                   # Logs de aplicacion, cron y procesos auxiliares
```

Uso de cada ubicacion:

- `/opt/containers_apps/<despliegue>/` agrupa una aplicacion o un conjunto de aplicaciones que se despliegan juntas. Contiene los clones Git, el codigo fuente y los ficheros de configuracion usados por la instalacion. La instalacion se ejecuta desde este directorio. Puede incluir un directorio `backup/` para backups propios del despliegue.
- `/srv/containers/backup/` es la ubicacion alternativa para centralizar los distintos tipos de backup. Cada despliegue debe elegir de forma coherente entre esta ruta y su directorio `backup/` bajo `/opt/containers_apps/`.
- `/srv/containers/bind/<aplicacion>/` contiene exclusivamente las rutas del host que se montan como bind mounts. Deben estar separadas por aplicacion y creadas con los propietarios y permisos requeridos antes de arrancar los contenedores. Estos permisos se gestionan mediante el script de instalación.
- `/srv/containers/storage/<usuario-podman>/` contiene la estructura de almacenamiento rootless de Podman, incluidos sus metadatos, capas, imagenes y volumenes. Podman gestiona esta estructura; no se deben cambiar manualmente sus propietarios o permisos. En los servidores actuales, `<usuario-podman>` puede ser, por ejemplo, `bioinfo`.
- `/srv/containers/shared/` se reserva para datos que deban compartir varias aplicaciones.
- `/var/log/local/<aplicacion>/` centraliza los logs persistentes del host. Como norma general, `apache/` contiene los logs del servidor web y `apps/` los de la aplicacion y sus procesos auxiliares. Algunas aplicaciones pueden necesitar subdirectorios adicionales.

Otros despliegues, como `beacon`, `localega` o `relecov-platform`, repiten esta misma separacion bajo su propio nombre. Esta convencion es importante al aplicar permisos: el codigo, los bind mounts, los volumenes gestionados por Podman y los logs no deben tratarse como si fueran el mismo tipo de almacenamiento.

## Preparar directorios del host

Los bind mounts son rutas reales del host. Deben existir y pertenecer al usuario que ejecuta `container_install.sh`.

Ejemplo recomendado:

```bash
sudo mkdir -p /var/log/local/relecov-iskylims/apps
sudo mkdir -p /var/log/local/relecov-iskylims/apache
sudo mkdir -p /srv/containers/bind/iskylims/iskylims_apache_conf
sudo mkdir -p /srv/containers/bind/iskylims/iskylims_django_setting

sudo chown -R "_USER-RUNNING_PODMAN_:_USER-RUNNING_PODMAN_" /var/log/local/relecov-iskylims
sudo chown -R "_USER-RUNNING_PODMAN_:_USER-RUNNING_PODMAN_" /srv/containers/bind/iskylims
```

Si usas otras rutas para `APACHE_CONF_PATH` o `DJANGO_SETTINGS_PATH`, crea esas rutas y asignales la misma propiedad.

`container_install.sh` ajustara despues los permisos internos con `podman unshare` y con `podman exec --user 0` cuando el contenedor este levantado.

## Actualizar codigo

Entra en el repositorio como el usuario que ejecuta Podman:

```bash
cd /ruta/al/repositorio/iskylims
git pull
```

Si todavia no existe el repositorio en el servidor:

```bash
git clone https://gitlab.isciii.es/bu-isciii/iSkyLIMS.git iskylims
cd iskylims
```

## Configurar `my_prod_settings_iskylims.txt`

Si el fichero ya existe, revisalo y mantenlo. Si no existe, copialo desde la plantilla:

```bash
cp conf/docker_production_settings.txt conf/my_prod_settings_iskylims.txt
```

Edita:

```bash
nano conf/my_prod_settings_iskylims.txt
```

Valores principales:

```bash
INSTALL_PATH='/opt/iskylims'

APACHE_CONF_PATH='/srv/containers/bind/iskylims/iskylims_apache_conf'
DJANGO_SETTINGS_PATH='/srv/containers/bind/iskylims/iskylims_django_setting/settings.py'

SERVER_STATUS_SERVER_NAME='<dns_server_status>'
SERVER_STATUS_ALIASES='127.0.0.1 localhost'
SERVER_STATUS_ALLOW_FROM='127.0.0.1 localhost'

APACHE_FORWARDED_PROTO='https'
APACHE_FORWARDED_PORT='443'

APP_UID='1212'
APP_GID='1212'
APP_SHELL='/sbin/nologin'
APP_PORT='8001'

DB_USER='<usuario_db>'
DB_PASS='<password_db>'
DB_NAME='iskylims'
DB_SERVER_IP='<host_o_ip_mysql>'
DB_PORT=3306

EMAIL_HOST_SERVER='docker.container.internal'
EMAIL_PORT='25'
EMAIL_HOST_USER='<correo>'
EMAIL_HOST_PASSWORD=''
EMAIL_USE_TLS='False'

LOCAL_SERVER_IP='*'
DNS_URL='<dns_o_ip_de_iskylims>'
```

Notas:

- `INSTALL_PATH` es la ruta dentro del contenedor.
- `APACHE_CONF_PATH` es una ruta del host donde se escriben los ficheros Apache renderizados.
- `DJANGO_SETTINGS_PATH` es una ruta del host para el `settings.py` montado en el contenedor.
- `DB_*` debe apuntar a la base de datos de produccion existente.
- Mantener el mismo `APP_UID` y `APP_GID` en todas las actualizaciones evita problemas de permisos.

## Backup antes de actualizar

Haz siempre backup antes de ejecutar la actualizacion.

Crea una carpeta de backup:

```bash
BACKUP_DIR=~/iskylims_backup_$(date +%Y%m%d_%H%M%S)
mkdir -p "$BACKUP_DIR"
```

Backup de base de datos:

```bash
mysqldump --user=<usuario_db> --password --host=<host_db> --port=<puerto_db> iskylims > "$BACKUP_DIR/iskylims.sql"
```

Backup de volumenes Podman:

```bash
podman volume ls | grep iskylims
podman volume export iskylims_iskylims_documents > "$BACKUP_DIR/iskylims_documents.tar"
podman volume export iskylims_iskylims_static > "$BACKUP_DIR/iskylims_static.tar"
```

Backup de configuracion:

```bash
cp conf/my_prod_settings_iskylims.txt "$BACKUP_DIR/"
cp .env.prod.file "$BACKUP_DIR/" 2>/dev/null || true
```

## Ejecutar la actualizacion

Para la mayoria de actualizaciones:

```bash
bash container_install.sh --engine podman --install_conf conf/my_prod_settings_iskylims.txt --action upgrade 2>&1 | tee ./iskylims_podman_upgrade_$(date +%Y%m%d_%H%M%S).log
```

El script:

- construye una nueva imagen;
- arranca o recrea los contenedores necesarios;
- genera `.env.prod.file`;
- renderiza configuracion Apache;
- prepara `settings.py`;
- repara permisos de bind mounts y volumenes;
- ejecuta `install.sh --bootstrap upgrade`;
- aplica migraciones;
- refresca `collectstatic`.

No usa la accion `install` porque este procedimiento asume una base de datos institucional ya existente.

## Caso especial: actualizacion desde 3.0.0 a 3.1.0

La actualizacion desde 3.0.0 a 3.1.0 requiere pasos extra porque hay cambios de datos que necesitan scripts especificos:

- `convert_rawtop_counter_to_int` antes de migraciones;
- `library_pool_to_many_relation` despues de migraciones;
- exportar antes la relacion `wetlab_library_pool.id -> run_process_id_id`.

Exporta el fichero necesario:

```bash
mysql --user=<usuario_db> --password --host=<host_db> --port=<puerto_db> iskylims \
  -e "SELECT id, run_process_id_id FROM wetlab_library_pool" \
  > /tmp/library_pool_run_process.tsv
```

Ejecuta la actualizacion especial:

```bash
bash container_install.sh --engine podman --install_conf conf/my_prod_settings_iskylims.txt --action upgrade \
  --script_before convert_rawtop_counter_to_int \
  --script_after library_pool_to_many_relation,/tmp/library_pool_run_process.tsv \
  2>&1 | tee ./iskylims_podman_upgrade_3_0_0_to_3_1_0_$(date +%Y%m%d_%H%M%S).log
```

Usa este comando solo para esa actualizacion concreta. Para actualizaciones posteriores, usa el comando generico de la seccion anterior.

## Comprobaciones posteriores

Comprueba contenedores:

```bash
podman compose --env-file .env.prod.file -f docker-compose.prod.yml ps
```

Revisa logs:

```bash
podman compose --env-file .env.prod.file -f docker-compose.prod.yml logs --tail 200 app
podman compose --env-file .env.prod.file -f docker-compose.prod.yml logs --tail 200 apache
```

Comprueba la aplicacion:

```text
http://<servidor>:8080
```

Si se han cambiado parametros de runtime en `conf/my_prod_settings_iskylims.txt`, vuelve a ejecutar `container_install.sh` para regenerar `.env.prod.file` y recrear los contenedores de forma coherente.

## Rollback

Si la actualizacion falla y necesitas volver atras:

1. Deten contenedores:

    ```bash
    podman compose --env-file .env.prod.file -f docker-compose.prod.yml down
    ```

2. Vuelve al commit o tag anterior del codigo:

    ```bash
    git checkout <commit_o_tag_anterior>
    ```

3. Restaura la base de datos:

    ```bash
    mysql --user=<usuario_db> --password --host=<host_db> --port=<puerto_db> iskylims < "$BACKUP_DIR/iskylims.sql"
    ```

4. Restaura volumenes si es necesario:

    ```bash
    podman volume import iskylims_iskylims_documents "$BACKUP_DIR/iskylims_documents.tar"
    podman volume import iskylims_iskylims_static "$BACKUP_DIR/iskylims_static.tar"
    ```

5. Restaura configuracion si cambio:

    ```bash
    cp "$BACKUP_DIR/my_prod_settings_iskylims.txt" conf/my_prod_settings_iskylims.txt
    ```

6. Repara permisos, arranca y vuelve a reparar volumenes montados:

    ```bash
    bash container_install.sh --engine podman --install_conf conf/my_prod_settings_iskylims.txt --action fix-permissions
    podman compose --env-file .env.prod.file -f docker-compose.prod.yml up -d
    bash container_install.sh --engine podman --install_conf conf/my_prod_settings_iskylims.txt --action fix-permissions
    ```

7. Revisa logs:

    ```bash
    podman compose --env-file .env.prod.file -f docker-compose.prod.yml logs --tail 200 app
    podman compose --env-file .env.prod.file -f docker-compose.prod.yml logs --tail 200 apache
    ```

## Reparar permisos

Ejecuta esta accion si:

- se han recreado contenedores manualmente;
- se han restaurado volumenes;
- se han cambiado propietarios en el host;
- se han cambiado `APP_UID` o `APP_GID`;
- el contenedor no arranca por permisos de bind mounts.

```bash
bash container_install.sh --engine podman --install_conf conf/my_prod_settings_iskylims.txt --action fix-permissions
```

Si el contenedor no esta arrancado, esta accion repara solo los bind mounts del host. Despues arranca los contenedores y repite la accion para reparar los volumenes montados:

```bash
bash container_install.sh --engine podman --install_conf conf/my_prod_settings_iskylims.txt --action fix-permissions
podman compose --env-file .env.prod.file -f docker-compose.prod.yml up -d
bash container_install.sh --engine podman --install_conf conf/my_prod_settings_iskylims.txt --action fix-permissions
```

## Operaciones utiles

Usa siempre `.env.prod.file` al ejecutar Podman Compose directamente:

```bash
podman compose --env-file .env.prod.file -f docker-compose.prod.yml ps
podman compose --env-file .env.prod.file -f docker-compose.prod.yml up -d
podman compose --env-file .env.prod.file -f docker-compose.prod.yml restart app
podman compose --env-file .env.prod.file -f docker-compose.prod.yml down
```

Entrar al contenedor:

```bash
podman exec -it iskylims_app bash
```

Ejecutar `collectstatic` manualmente:

```bash
podman exec -it iskylims_app bash -lc 'cd /opt/iskylims && source virtualenv/bin/activate && python manage.py collectstatic --noinput'
```

Ejecutar manualmente el bootstrap de actualizacion:

```bash
podman exec -it iskylims_app bash -c 'cd /srv/iskylims && bash install.sh --bootstrap upgrade --git_revision main --conf conf/my_prod_settings_iskylims.txt --tables --skip_apache_restart'
```

## Notas de permisos

Bind mounts:

- Son rutas reales del host.
- Deben existir antes de arrancar contenedores.
- Deben pertenecer al usuario que ejecuta Podman rootless.
- `container_install.sh` usa `podman unshare` para aplicar propietarios internos cuando hace falta.

Volumenes Podman:

- Los gestiona Podman en el almacenamiento rootless del usuario.
- Se reparan desde dentro del contenedor con `podman exec --user 0`.
- Si cambias `APP_UID` o `APP_GID`, ejecuta `--action fix-permissions`.

Apache:

- El contenedor Apache UBI usa UID `1001` y grupo `0`.
- Los logs Apache se preparan para ese usuario.
- Los ficheros Apache renderizados se dejan con permisos `0664`.

`settings.py`:

- Se monta desde el host.
- `container_install.sh` lo prepara con permisos `0664`.
- Si editas el fichero a mano, ejecuta despues `--action fix-permissions`.
