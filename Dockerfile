FROM ubuntu:24.04
ENV TZ=Europe/Madrid
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Updates
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get upgrade -y

# Essential software
RUN apt-get install -y \
    git wget lsb-release \
    libmysqlclient-dev default-mysql-client \
    python3-pip libpq-dev python3-venv python3-wheel \
    apache2-dev cron \
    gnuplot pkg-config rsync

# Set MYSQLCLIENT_CFLAGS and MYSQLCLIENT_LDFLAGS using pkg-config
RUN export MYSQLCLIENT_CFLAGS="$(pkg-config --libs mysqlclient)" && \
    export MYSQLCLIENT_LDFLAGS="$(pkg-config --cflags mysqlclient)"

# Set git repository
RUN mkdir /srv/iskylims 
WORKDIR /srv/iskylims

# Copy the local git repository to docker image directory
COPY . /srv/iskylims

ENV PATH="/usr/sbin/cron:$PATH"

# Set default install type
ARG INSTALL_TYPE=dep
ARG GIT_REVISION=main
ARG INSTALL_CONF=conf/docker_install_settings.txt

# Execute the dependency stage only; app migrations run when the container is up.
RUN /bin/bash install.sh --install dep --git_revision $GIT_REVISION --conf $INSTALL_CONF --skip_apache_restart
# Use the virtualenv created by install.sh
ENV PATH="/opt/relecov-platform/virtualenv/bin:${PATH}"

WORKDIR /opt/iskylims

# Expose
EXPOSE 8001

# Start the application
CMD ["python", "/opt/iskylims/manage.py", "runserver", "0:8001"]
