FROM registry.access.redhat.com/ubi9/ubi
ENV TZ=Europe/Madrid
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Runtime user (override with build args if needed)
ARG APP_UID=1212
ARG APP_GID=1212


# Updates
RUN dnf -y update

# Essential software
RUN dnf -y install \
    git wget \
    python3.11 python3.11-pip python3.11-devel python3.11-wheel \
    gcc gcc-c++ make \
    openssl-devel libffi-devel \
    mariadb mariadb-connector-c-devel postgresql-devel \
    httpd-devel cronie \
    rsync tzdata \
    pkgconf-pkg-config \
    && dnf clean all

# Ensure python3 points to the desired version
RUN ln -sf /usr/bin/python3.11 /usr/bin/python3

# Set MYSQLCLIENT_CFLAGS and MYSQLCLIENT_LDFLAGS using pkg-config
RUN export MYSQLCLIENT_CFLAGS="$(pkg-config --libs mysqlclient)" && \
    export MYSQLCLIENT_LDFLAGS="$(pkg-config --cflags mysqlclient)"

# Set git repository
RUN mkdir /srv/iskylims 
WORKDIR /srv/iskylims

# Copy the local git repository to docker image directory
COPY . /srv/iskylims

ENV PATH="/usr/sbin/cron:$PATH"
RUN chmod +x /srv/iskylims/scripts/container_start.sh

# Set default install type
ARG INSTALL_TYPE=dep
ARG GIT_REVISION=main
ARG INSTALL_CONF=conf/docker_test_settings.txt

# Execute the dependency stage only; app migrations run when the container is up.
ENV SKIP_SYSTEM_PACKAGES=1
RUN /bin/bash install.sh --install dep --git_revision $GIT_REVISION --conf $INSTALL_CONF --skip_apache_restart
# Use the virtualenv created by install.sh
ENV PATH="/opt/iskylims/virtualenv/bin:${PATH}"

WORKDIR /opt/iskylims

# Create non-root user and set ownership
RUN groupadd -g ${APP_GID} iskylims && \
    useradd -m -u ${APP_UID} -g ${APP_GID} -s /sbin/nologin iskylims && \
    chown -R ${APP_UID}:${APP_GID} /opt/iskylims /srv/iskylims && \
    git config --system --add safe.directory /srv/iskylims

# Expose
EXPOSE 8001

# Start the application once install.sh has populated /opt/iskylims.
USER iskylims
CMD ["/srv/iskylims/scripts/container_start.sh"]
