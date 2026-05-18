FROM registry.access.redhat.com/ubi9/ubi
ENV TZ=Europe/Madrid
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Runtime user (override with build args if needed)
ARG APP_UID=1212
ARG APP_GID=1212
ARG APP_SHELL=/sbin/nologin
ARG APP_INSTALL_PATH=/opt/iskylims
ENV APP_INSTALL_PATH=${APP_INSTALL_PATH}
ENV PIP_NO_CACHE_DIR=1

# Updates
RUN dnf -y update \
    && dnf clean all \
    && rm -rf /var/cache/dnf /tmp/* /var/tmp/*

# Add EPEL for packages not available in default UBI repositories
RUN dnf -y install --setopt=install_weak_deps=False --nodocs \
    https://dl.fedoraproject.org/pub/epel/epel-release-latest-9.noarch.rpm \
    && dnf clean all \
    && rm -rf /var/cache/dnf /tmp/* /var/tmp/*

# Essential software
RUN dnf -y install --setopt=install_weak_deps=False --nodocs \
    git wget \
    python3.11 python3.11-pip python3.11-devel python3.11-wheel \
    gcc gcc-c++ make \
    openssl-devel libffi-devel \
    mariadb mariadb-connector-c-devel \
    httpd-devel cronie \
    rsync tzdata \
    pkgconf-pkg-config \
    gnuplot-minimal \
    && dnf clean all \
    && rm -rf /var/cache/dnf /tmp/* /var/tmp/*

# Install supercronic (rootless-friendly cron runner)
RUN set -eux; \
    SUPERCRONIC_VERSION="v0.2.38"; \
    arch="$(uname -m)"; \
    case "$arch" in \
      x86_64) supercronic_arch="amd64" ;; \
      aarch64) supercronic_arch="arm64" ;; \
      *) echo "Unsupported architecture for supercronic: $arch" >&2; exit 1 ;; \
    esac; \
    supercronic_url="https://github.com/aptible/supercronic/releases/download/${SUPERCRONIC_VERSION}/supercronic-linux-${supercronic_arch}"; \
    if wget --tries=3 --waitretry=2 --retry-connrefused -q -O /usr/local/bin/supercronic "${supercronic_url}"; then \
      chmod +x /usr/local/bin/supercronic; \
    else \
      rm -f /usr/local/bin/supercronic; \
      echo "supercronic download failed from ${supercronic_url}; continuing without cron support"; \
    fi; \
    rm -rf /tmp/* /var/tmp/*

# Ensure python3 points to the desired version
RUN ln -sf /usr/bin/python3.11 /usr/bin/python3

# Install Illumina InterOp CLI used to generate run metric plots
RUN set -eux; \
    cd /opt; \
    wget -q https://github.com/Illumina/interop/releases/download/v1.1.15/InterOp-1.1.15-Linux-GNU.tar.gz; \
    tar -xf InterOp-1.1.15-Linux-GNU.tar.gz; \
    ln -s InterOp-1.1.15-Linux-GNU interop; \
    rm InterOp-1.1.15-Linux-GNU.tar.gz; \
    rm -rf /tmp/* /var/tmp/*

# Set git repository
RUN mkdir /srv/iskylims 
WORKDIR /srv/iskylims

# Copy the local git repository to docker image directory
COPY --chown=${APP_UID}:${APP_GID} . /srv/iskylims

ENV PATH="/usr/sbin/cron:$PATH"
RUN chmod +x /srv/iskylims/scripts/container_start.sh

# Set default install type
ARG INSTALL_TYPE=dep
ARG GIT_REVISION=main
ARG INSTALL_CONF=conf/docker_test_settings.txt

# Prepare dependencies and stage the application tree in the image so the
# container can restart without rerunning install-time file generation.
ENV SKIP_SYSTEM_PACKAGES=1
RUN /bin/bash install.sh --install dep --git_revision $GIT_REVISION --conf $INSTALL_CONF --skip_apache_restart \
    && rm -rf /root/.cache/pip /tmp/* /var/tmp/*
RUN /bin/bash install.sh --stage install --git_revision $GIT_REVISION --conf $INSTALL_CONF --skip_apache_restart \
    && rm -rf /root/.cache/pip /tmp/* /var/tmp/*
# Use the virtualenv created by install.sh
ENV PATH="${APP_INSTALL_PATH}/virtualenv/bin:${PATH}"

WORKDIR ${APP_INSTALL_PATH}

# Create non-root user and set ownership
RUN groupadd -g ${APP_GID} iskylims && \
    useradd -m -u ${APP_UID} -g ${APP_GID} -s ${APP_SHELL} iskylims && \
    mkdir -p ${APP_INSTALL_PATH}/cron ${APP_INSTALL_PATH}/tmp && \
    chown -R ${APP_UID}:${APP_GID} ${APP_INSTALL_PATH} /srv/iskylims && \
    chmod 700 ${APP_INSTALL_PATH}/cron ${APP_INSTALL_PATH}/tmp && \
    git config --system --add safe.directory /srv/iskylims

# Expose
EXPOSE 8001

# Start the application once install.sh has populated /opt/iskylims.
USER iskylims
CMD ["/srv/iskylims/scripts/container_start.sh"]
