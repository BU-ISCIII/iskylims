# iSkyLIMS

[![Django](https://img.shields.io/static/v1?label=Django&message=4.2&color=azul?style=plastic&logo=django)](https://github.com/django/django)
[![Python](https://img.shields.io/static/v1?label=Python&message=3.8.10&color=verde?style=plastic&logo=Python)](https://www.python.org/)
[![Bootstrap](https://img.shields.io/badge/Bootstrap-v5.0-azulvioleta?style=plastic&logo=Bootstrap)](https://getbootstrap.com)
[![version](https://img.shields.io/badge/version-3.0.0-naranja?style=plastic&logo=GitHub)](https://github.com/BU-ISCIII/iskylims.git)

La introduccion de la secuenciacion masiva (MS) en las instalaciones de genomica ha significado un crecimiento exponencial en la generacion de datos, lo que requiere un sistema de seguimiento preciso, desde la preparacion de la biblioteca hasta la generacion de archivos fastq, el analisis y la entrega al investigador. El software disenado para manejar esas tareas se llama Sistemas de Gestion de Informacion de Laboratorio (LIMS), y su software debe adaptarse a las necesidades particulares de su laboratorio de genomica. iSkyLIMS nace con el objetivo de ayudar con las tareas de laboratorio humedo e implementar un flujo de trabajo que guie a los laboratorios de genomica en sus actividades, desde la preparacion de la biblioteca hasta la produccion de datos, reduciendo los posibles errores asociados a la tecnologia de alto rendimiento y facilitando el control de calidad de la secuenciacion. Ademas, iSkyLIMS conecta el laboratorio humedo con el laboratorio seco, facilitando el analisis de datos por parte de bioinformaticos.

![Imagen](img/iskylims_scheme.png)

De acuerdo con la infraestructura existente, la secuenciacion se realiza en un instrumento Illumina NextSeq. Los datos se almacenan en un dispositivo de almacenamiento masivo NetApp y los archivos fastq se generan (bcl2fastq) en un cluster de computo de alto rendimiento Sun Grid Engine (SGE-HPC). Los servidores de aplicaciones ejecutan aplicaciones web para el analisis bioinformatico (GALAXY), la aplicacion iSkyLIMS y alojan la capa de informacion de MySQL. El flujo de trabajo de iSkyLIMS WetLab se ocupa del seguimiento y las estadisticas de la ejecucion de la secuenciacion. El seguimiento de la ejecucion pasa por cinco estados: "registrado", el usuario de genomica registra la nueva ejecucion de la secuenciacion en el sistema, el proceso esperara hasta que la ejecucion se complete en la maquina y los datos se transfieran al dispositivo de almacenamiento masivo; "Envio de hoja de muestra", el archivo de hoja de muestra con la informacion de la ejecucion de la secuenciacion se copiara en la carpeta de ejecucion para el proceso de bcl2fastq; "Procesamiento de datos", se procesan los archivos de parametros de ejecucion y los datos se almacenan en la base de datos; "Estadisticas en ejecucion", los datos de desmultiplexacion generados en el proceso de bcl2fastq se procesan y almacenan en la base de datos, "Completado", todos los datos se procesan y almacenan correctamente. Se proporcionan estadisticas por muestra, por proyecto, por ejecucion y por investigacion, asi como informes anuales y mensuales. El flujo de trabajo de iSkyLIMS DryLab se encarga de la solicitud de servicios de bioinformatica y estadisticas. El usuario solicita servicios que pueden estar asociados con una ejecucion de secuenciacion. Se proporciona seguimiento de estadisticas y servicios.

- [iSkyLIMS](#iskylims)
  - [Obtener el codigo (obligatorio)](#obtener-el-codigo-obligatorio)
  - [Elige tu ruta](#elige-tu-ruta)
  - [Despliegue con Docker](#despliegue-con-docker)
    - [Contenedor local de pruebas](#contenedor-local-de-pruebas)
    - [Contenedor de produccion](#contenedor-de-produccion)
      - [Proxy inverso con Apache (host) + Gunicorn](#proxy-inverso-con-apache-host--gunicorn)
    - [Actualizacion del despliegue Docker](#actualizacion-del-despliegue-docker)
    - [Actualizacion del despliegue Docker v3.0.0 a 3.1.0](#actualizacion-del-despliegue-docker-v300-a-310)
      - [Haz copia de seguridad](#haz-copia-de-seguridad)
  - [Despliegue bare-metal (Ubuntu/CentOS)](#despliegue-bare-metal-ubuntucentos)
    - [Instalacion](#instalacion)
      - [Requisitos previos](#requisitos-previos)
      - [Clonar el repositorio](#clonar-el-repositorio)
      - [Preparar la base de datos](#preparar-la-base-de-datos)
      - [Configurar install\_settings.txt](#configurar-install_settingstxt)
      - [Ejecutar install.sh](#ejecutar-installsh)
    - [Actualizacion (3.0.x a 3.1.x)](#actualizacion-30x-a-31x)
      - [Haz copia de seguridad](#haz-copia-de-seguridad-1)
      - [Actualizar codigo y ajustes](#actualizar-codigo-y-ajustes)
      - [Ejecutar pasos de actualizacion con root](#ejecutar-pasos-de-actualizacion-con-root)
      - [Ejecutar pasos de actualizacion sin root](#ejecutar-pasos-de-actualizacion-sin-root)
  - [Que hacer si algo falla](#que-hacer-si-algo-falla)
    - [Bare-metal](#bare-metal)
    - [Docker](#docker)
  - [Pasos finales de configuracion](#pasos-finales-de-configuracion)
    - [Configuracion de SAMBA](#configuracion-de-samba)
    - [Verificacion de correo electronico](#verificacion-de-correo-electronico)
    - [Configurar el servidor Apache](#configurar-el-servidor-apache)
    - [Verificacion de la instalacion](#verificacion-de-la-instalacion)
  - [Documentacion de iSkyLIMS](#documentacion-de-iskylims)

Si tienes algun problema o deseas informar de algun error, por favor, publicalo en [issue](https://github.com/BU-ISCIII/iSkyLIMS/issues)

## Obtener el codigo (obligatorio)

Todas las rutas de instalacion asumen que ya clonaste el repositorio:

```bash
git clone https://github.com/BU-ISCIII/iskylims.git iskylims
cd iskylims
```

## Elige tu ruta

- **Docker (pruebas locales)**: levanta MySQL + Samba + iSkyLIMS con datos de demo para probar rapidamente.
- **Docker (contenedor de produccion)**: despliega solo la aplicacion, apuntando a tu DB/Samba existente.
- **Bare-metal**: instala o actualiza directamente en hosts Ubuntu/CentOS con `install.sh`.

## Despliegue con Docker

Requisitos previos para instalaciones con Docker:

- Docker Engine + Docker Compose v2
- git (para clonar el repositorio)

### Contenedor local de pruebas

Levanta el sistema completo (base de datos, Samba y app) con fixtures y datos de demo:

```bash
bash docker_install.sh --test
```

Esto usa `docker-compose.test.yml` por defecto.

Puedes personalizar los valores por defecto:

- `--demo_data /ruta/a/iskylims_demo_data.tar.gz` para reutilizar un archivo local (si no, se descarga).
- `--skip_demo_data` o `--skip_test_data` para evitar cargar datos extra.
- `--install_type` (`full` por defecto) y `--git_revision` para controlar el build.
- `--script` para ejecutar uno o mas scripts de migracion via `install.sh` (puedes repetir la opcion).

Ejemplo de uso con un script de migracion en Docker:

```bash
bash docker_install.sh --test --script migrate_optional_values
```

Cuando el script termine, abre `http://localhost:8001` y crea el superusuario de Django cuando te lo pida.

### Contenedor de produccion

Despliega el contenedor de iSkyLIMS contra servicios MySQL/Samba externos:

1. Copia y edita la plantilla de produccion:

    ```bash
    cp conf/docker_production_settings.txt conf/my_prod_settings.txt
    # edita conf/my_prod_settings.txt con tus datos de DB/Samba
    ```

2. Construye y ejecuta en modo produccion (usa `docker-compose.prod.yml` por defecto):

    ```bash
    bash docker_install.sh --install_conf conf/my_prod_settings.txt
    ```

   Usa `--compose_file` para cambiar el compose o `--install_type`/`--git_revision` para variar el build.
   Tip: captura logs para depuracion:

    ```bash
    bash docker_install.sh --install_conf conf/my_prod_settings.txt 2>&1 | tee ./iskylims_docker_install_$(date +%Y%m%d_%H%M%S).log
    ```

3. Si es una instalacion nueva, crea el superusuario cuando se solicite y completa la configuracion de Samba en la UI.

UID/GID del usuario de ejecucion del contenedor (por defecto `1212:1212`):

- Exporta `APP_UID` y `APP_GID` antes de ejecutar `docker_install.sh` si necesitas un UID/GID distinto (por ejemplo, para escribir en `/opt/iskylims/static-host`).

```bash
export APP_UID=1212
export APP_GID=1212
```

Asegura que la carpeta de estaticos en el host sea escribible por ese UID/GID:

```bash
sudo chown -R 1212:1212 /opt/iskylims/static-host
```

#### Proxy inverso con Apache (host) + Gunicorn

En produccion, el contenedor ejecuta `gunicorn` (no `manage.py runserver`). Usa Apache en el host como proxy inverso hacia `localhost:8001`.

#### Tareas cron dentro del contenedor

Cron se ejecuta mediante un `crond` ligero lanzado por el script de arranque del contenedor. El script escribe las entradas de django-crontab en `/opt/iskylims/cron/iskylims` y arranca `crond` con un PID en una ruta escribible por el usuario.

Si modificas `CRONJOBS`, reconstruye o reinicia el contenedor para regenerar el archivo de cron.

Archivos estaticos:

- El contenedor genera los estaticos en `/opt/iskylims/static`.
- `docker-compose.prod.yml` monta ese directorio en `/opt/iskylims/static-host` del host.
- Configura Apache con `Alias /static/ /opt/iskylims/static-host/`.

Ejemplo de configuracion en `conf/iskylims_apache_reverse_proxy.conf`.

### Actualizacion del despliegue Docker

Re-despliega el contenedor de aplicacion contra una base de datos existente sin tocar los datos:

```bash
bash docker_install.sh --install_conf conf/my_prod_settings.txt --action upgrade
```

La actualizacion reconstruye/reinicia el contenedor y ejecuta `install.sh` dentro del contenedor, que regenera migraciones, las aplica con `--fake-initial` y evita cargar superusuario/datos demo/prueba.

Si usas `APP_UID`/`APP_GID`, exportalos de nuevo antes de la actualizacion para mantener el mismo UID/GID:

```bash
export APP_UID=1212
export APP_GID=1212
```

### Actualizacion del despliegue Docker v3.0.0 a 3.1.0

#### Haz copia de seguridad

- Copia completa de la base de datos `iskylims`.
- Copia completa de las carpetas de logs y documents.

```bash
docker run --rm \
  -v iskylims_logs:/from \
  -v "$PWD":/to \
  alpine tar -czf /to/iskylims_logs.tgz -C /from .

docker run --rm \
  -v iskylims_documents:/from \
  -v "$PWD":/to \
  alpine tar -czf /to/iskylims_documents.tgz -C /from .
```

Antes de actualizar, asegurate de tener una copia completa de la base de datos y de los datos en volumenes que uses. Confirma que `conf/my_prod_settings.txt` tenga el host/usuario/password de la base de datos de produccion, la URL/IP del servidor, correo y ajustes de logging usados por el contenedor.

Para 3.0.0 -> 3.1.0, exporta primero el mapeo de LibraryPool y luego ejecuta la actualizacion con scripts pre/post:

```bash
mysql --user=<db_user> --password=<db_password> --host=<db_server_ip> --port=<db_port> iskylims \
  -e "SELECT id, run_process_id_id FROM wetlab_library_pool" \
  > /tmp/library_pool_run_process.tsv

bash docker_install.sh --install_conf conf/my_prod_settings.txt --action upgrade \
  --script_before convert_rawtop_counter_to_int \
  --script_after library_pool_to_many_relation,/tmp/library_pool_run_process.tsv
```

## Despliegue bare-metal (Ubuntu/CentOS)

### Instalacion

#### Requisitos previos

- **Privilegios sudo** para instalar dependencias
- MySQL > 8.0 o MariaDB > 10.4
- Apache 2.4
- git > 2.34
- Python > 3.11
- Servidor local configurado para enviar correos
- Acceso a la carpeta Samba donde estan los run folders
- Paquete `lsb_release` (`yum install redhat-lsb-core` en RedHat/CentOS, `apt install lsb-core lsb-release` en Ubuntu)

#### Clonar el repositorio

```bash
cd <tu directorio de trabajo>
git clone https://github.com/BU-ISCIII/iskylims.git iskylims
cd iskylims
```

#### Preparar la base de datos

1. Crea una base de datos llamada `iskylims`.
2. Crea un usuario con permisos de lectura/escritura sobre esa base.
3. Guarda host, puerto, usuario y password para el archivo de ajustes.

#### Configurar install_settings.txt

```bash
cp conf/template_install_settings.txt install_settings.txt
nano install_settings.txt
```

Completa los valores de base de datos, email, IP/URL del servidor y logging.

#### Ejecutar install.sh

iSkyLIMS se instala en `/opt/iskylims` por defecto. El script `install.sh` gestiona dependencias y aplicacion; elige lo que necesitas con `--install`:

- `dep`: instala dependencias del sistema y de Python (requiere sudo).
- `app`: despliega el codigo, actualiza ajustes, ejecuta migraciones y collectstatic (sin sudo).
- `full`: ejecuta ambos pasos.

Ejemplos:

```bash
# solo dependencias del sistema
sudo bash install.sh --install dep

# solo aplicacion iSkyLIMS
bash install.sh --install app --git_revision main --tables

# dependencias + aplicacion
sudo bash install.sh --install full --git_revision main --tables
```

- Añade `--tables` para cargar los datos iniciales en instalaciones nuevas, o `--skip_tables` si quieres omitirlos.
- Captura logs para depuracion con `tee`:

  ```bash
  sudo bash install.sh --install full --git_revision main --tables 2>&1 | tee ./iskylims_install_$(date +%Y%m%d_%H%M%S).log
  ```

- Si Apache se gestiona desde otro sitio, omite el reinicio automatico con `--skip_apache_restart`.

### Actualizacion (3.0.x a 3.1.x)

Sigue estos pasos para pasar de la version 3.0.0 a la serie 3.1.x.

#### Haz copia de seguridad

- Copia completa de la base de datos `iskylims`.
- Copia completa de la carpeta de instalacion (por ejemplo `/opt/iskylims`).
- Si usas library pools, exportalos antes de actualizar:

  ```bash
  mysql --user=<db_user> --password=<db_password> --host=<db_server_ip> --port=<db_port> iskylims \
    -e "SELECT * FROM wetlab_library_pool" > <carpeta_backup>/backup_lib_pool.sql
  ```

#### Actualizar codigo y ajustes

```bash
cd <tu directorio de trabajo>/iskylims
git pull
cp conf/template_install_settings.txt install_settings.txt
sudo nano install_settings.txt
```

Si editas el archivo en Windows, asegurate de guardarlo con codificacion UTF-8/ASCII.

#### Ejecutar pasos de actualizacion con root

Actualiza dependencias del sistema y de Python:

```bash
sudo bash install.sh --upgrade dep 2>&1 | tee install_full.log
```

Asegura que los permisos permiten que el paso sin root escriba en `/opt/iskylims` (ajusta tu hardening si cambio la ruta).

#### Ejecutar pasos de actualizacion sin root

Actualiza el codigo y la base de datos:

```bash
# con restauracion de library pool
bash install.sh --upgrade app --script <carpeta_backup>/backup_lib_pool.sql --git_revision main --tables

# sin restauracion de library pool
bash install.sh --upgrade app --git_revision main --tables

# ejemplo ejecutando un script de migracion en la actualizacion
bash install.sh --upgrade app --script migrate_optional_values --git_revision main --tables

# 3.0.0 -> 3.1.0 (scripts de datos pre/post)
mysql --user=<db_user> --password=<db_password> --host=<db_server_ip> --port=<db_port> iskylims \
  -e "SELECT id, run_process_id_id FROM wetlab_library_pool" \
  > /tmp/library_pool_run_process.tsv

bash install.sh --upgrade app --git_revision main \
  --script_before convert_rawtop_counter_to_int \
  --script_after library_pool_to_many_relation,/tmp/library_pool_run_process.tsv
```

O ejecuta todo en un unico comando:

```bash
sudo bash install.sh --upgrade full --git_revision main --tables
```

Las actualizaciones regeneran las migraciones y las aplican con `--fake-initial` para conservar las tablas existentes, igual que en Docker.

## Que hacer si algo falla

Cuando actualizamos usando el script de instalacion estamos realizando varios cambios en la base de datos. Si algo falla necesitamos restaurar el estado anterior y empezar de nuevo.

### Bare-metal

Necesitamos copiar la carpeta completa `/opt/iskylims` de vuelta a `/opt/iskylims` (o tu ruta de instalacion), y restaurar la base de datos con algo como:

```bash
sudo rm -rf /opt/iskylims
sudo cp -r /home/dadmin/backup_prod/iSkyLIMS/ /opt/
sudo /scripts/hardening.sh
mysql -u iskylims -p -h dmysqlps.isciiides.es
# drop database iskylims;
# create database iskylims;
mysql -u iskylims -p -h dmysqlps.isciiides.es iskylims < /home/dadmin/backup_prod/bk_iSkyLIMS_202310160737.sql
```

### Docker

1. Para el contenedor:

    ```bash
    docker compose -f docker-compose.prod.yml down
    ```

2. Restaura la base de datos desde tu backup.

3. Si sospechas que la imagen o la cache de build esta corrupta, elimina la imagen de la app y reconstruye:

    ```bash
    docker image ls | grep iskylims
    docker rmi <iskylims_image_id>
    ```

4. Si los volumenes estan comprometidos, restauralos desde los tar:

    ```bash
    docker run --rm -v iskylims_logs:/to -v "$PWD":/from alpine \
      tar -xzf /from/iskylims_logs.tgz -C /to

    docker run --rm -v iskylims_documents:/to -v "$PWD":/from alpine \
      tar -xzf /from/iskylims_documents.tgz -C /to
    ```

5. Arranca el contenedor de nuevo:

    ```bash
    bash docker_install.sh --install_conf conf/my_prod_settings.txt --action upgrade
    ```

## Pasos finales de configuracion

### Configuracion de SAMBA

- Inicia sesion con la cuenta admin.
- Ve a Massive sequencing
![go_to_wetlab](img/got_to_wetlab.png){width:50px}
- Ve a Configuration -> Samba configuration
- Rellena el formulario con los parametros apropiados para la carpeta compartida de Samba:
![samba form](img/samba_form.png)

### Verificacion de correo electronico

- Ve a Massive sequencing
- Ve a Configuration -> Email configuration
- Rellena el formulario con los parametros necesarios y prueba a enviar un correo.

### Configurar el servidor Apache

Copia el archivo de configuracion de Apache segun tu distribucion dentro del directorio de configuracion de Apache y renombralo a iskylims.conf

Ubicaciones tipicas:

- Ubuntu/Debian: `/etc/apache2/sites-available/iskylims.conf` (habilitar con `a2ensite`)
- CentOS/RHEL: `/etc/httpd/conf.d/iskylims.conf`

Pasos sugeridos (Apache en el host como proxy inverso):

1. Copia el ejemplo de configuracion:

    ```bash
    sudo cp conf/iskylims_apache_reverse_proxy.conf /etc/apache2/sites-available/iskylims.conf
    # CentOS/RHEL:
    # sudo cp conf/iskylims_apache_reverse_proxy.conf /etc/httpd/conf.d/iskylims.conf
    ```

2. Edita la configuracion:

    - Ajusta `ServerName`
    - Comprueba que `ProxyPass` apunte a `http://localhost:8001/`
    - Comprueba `Alias /static/ /opt/iskylims/static-host/`

3. Crea la carpeta de estaticos en el host:

    ```bash
    sudo mkdir -p /opt/iskylims/static-host
    ```

4. Habilita modulos necesarios (Ubuntu/Debian):

    ```bash
    sudo a2enmod proxy proxy_http headers
    sudo a2ensite iskylims.conf
    ```

5. Recarga Apache:

    ```bash
    sudo systemctl reload apache2
    # CentOS/RHEL:
    # sudo systemctl reload httpd
    ```

### Verificacion de la instalacion

Abre el navegador y escribe "localhost" o la IP local del servidor para comprobar que iSkyLIMS esta funcionando.

Tambien puedes comprobar parte de la funcionalidad y las conexiones a Samba y base de datos usando:

- Ve a [configuration test](https://iskylims.isciii.es/wetlab/configurationTest/)
- Haz click en submit
- Revisa todas las pestañas para confirmar que la conexion es correcta.
- Ejecuta las 3 pruebas para cada maquina de secuenciacion: MiSeq, NextSeq y NovaSeq.

## Documentacion de iSkyLIMS

La documentacion de iSkyLIMS esta disponible en [https://iskylims.readthedocs.io/en/latest](https://iskylims.readthedocs.io/en/latest)
