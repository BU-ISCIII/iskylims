# syntax=docker/dockerfile:1.4
FROM registry.access.redhat.com/ubi9/ubi-minimal

ARG APP_UID=1212
ARG APP_GID=1212
ARG APP_PORT=8001
ARG APP_REPO_PATH=/srv/iskylims
ARG APP_INSTALL_PATH=/opt/iskylims
ARG GIT_REVISION=current
ARG INSTALL_CONF=conf/docker_test_settings.txt
ARG USE_INSTALL_CONF_SECRET=false
ARG RENDER_DJANGO_SETTINGS=false

ENV APP_REPO_PATH=${APP_REPO_PATH} \
    APP_INSTALL_PATH=${APP_INSTALL_PATH} \
    INSTALL_PATH=${APP_INSTALL_PATH} \
    APP_PORT=${APP_PORT} \
    PROJECT_MODULE=iskylims \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    TZ=Europe/Madrid

RUN microdnf -y update && \
    microdnf -y install python3.11 python3.11-pip \
      python3.11-devel tar gcc git rsync wget mariadb-connector-c-devel shadow-utils && \
    microdnf clean all

ARG SUPERCRONIC_VERSION=v0.2.38
RUN set -eux; \
    arch="$(uname -m)"; \
    case "$arch" in \
      x86_64) supercronic_arch="amd64" ;; \
      aarch64) supercronic_arch="arm64" ;; \
      *) echo "Unsupported architecture for supercronic: $arch" >&2; exit 1 ;; \
    esac; \
    wget -q -O /usr/local/bin/supercronic \
      "https://github.com/aptible/supercronic/releases/download/${SUPERCRONIC_VERSION}/supercronic-linux-${supercronic_arch}"; \
    chmod 0755 /usr/local/bin/supercronic

WORKDIR ${APP_REPO_PATH}
COPY . ${APP_REPO_PATH}/
COPY scripts/container_start.sh /usr/local/bin/container_start.sh

# Production uses an ephemeral build-secret mount: the operator configuration
# is readable by install.sh but never copied into an image layer. Test builds
# use the committed safe profile and explicitly render test settings.
RUN --mount=type=secret,id=install_conf \
    conf_path="${INSTALL_CONF}"; \
    if [ "${USE_INSTALL_CONF_SECRET}" = "true" ]; then \
      conf_path=/run/secrets/install_conf; \
      test -f "$conf_path" || { echo "Required install_conf build secret is missing" >&2; exit 1; }; \
    fi; \
    render_args=""; \
    if [ "${RENDER_DJANGO_SETTINGS}" = "true" ]; then render_args="--render-settings"; fi; \
    bash install.sh --stage install \
      --git_revision "${GIT_REVISION}" --conf "$conf_path" $render_args && \
    groupadd -g "${APP_GID}" app && \
    useradd -u "${APP_UID}" -g "${APP_GID}" -s /sbin/nologin app && \
    chmod 0755 /usr/local/bin/container_start.sh && \
    chown -R "${APP_UID}:${APP_GID}" "${APP_INSTALL_PATH}" "${APP_REPO_PATH}"

WORKDIR ${APP_INSTALL_PATH}
USER ${APP_UID}:${APP_GID}
EXPOSE ${APP_PORT}
HEALTHCHECK --interval=30s --timeout=5s --start-period=30s --retries=3 \
  CMD python -c "import os, urllib.request; urllib.request.urlopen('http://127.0.0.1:' + os.environ['APP_PORT'] + '/health/', timeout=3)" || exit 1
CMD ["/usr/local/bin/container_start.sh"]
