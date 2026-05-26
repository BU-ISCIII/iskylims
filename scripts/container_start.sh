#!/usr/bin/env bash
set -euo pipefail

APP_DIR="${INSTALL_PATH:-/opt/iskylims}"
CRON_DIR="${APP_DIR}/cron"
TMP_DIR="${APP_DIR}/tmp"
CRON_FILE="${CRON_DIR}/iskylims"
CRON_LOG="${TMP_DIR}/supercronic.log"
APP_MODE="${APP_MODE:-prod}"
APP_PORT="${APP_PORT:-8001}"
PROJECT_MODULE="${PROJECT_MODULE:-iskylims}"
GUNICORN_TIMEOUT="${GUNICORN_TIMEOUT:-120}"
GUNICORN_KEEPALIVE="${GUNICORN_KEEPALIVE:-5}"
GUNICORN_THREADS="${GUNICORN_THREADS:-2}"

if [ ! -f "${APP_DIR}/manage.py" ]; then
    echo "Application entrypoint not found at ${APP_DIR}/manage.py" >&2
    ls -la "${APP_DIR}" >&2 || true
    exit 1
fi

if [ ! -f "${APP_DIR}/virtualenv/bin/activate" ]; then
    echo "Virtualenv activation script not found at ${APP_DIR}/virtualenv/bin/activate" >&2
    ls -la "${APP_DIR}/virtualenv" >&2 || true
    exit 1
fi

source "${APP_DIR}/virtualenv/bin/activate"

mkdir -p "${CRON_DIR}" "${TMP_DIR}"

safe_chmod() {
    local mode="$1"
    shift
    local path
    for path in "$@"; do
        if [ ! -e "${path}" ]; then
            continue
        fi
        if [ -O "${path}" ]; then
            chmod "${mode}" "${path}"
        else
            echo "Skipping chmod ${mode} on ${path}: not owned by $(id -un)."
        fi
    done
}

safe_chmod 700 "${CRON_DIR}" "${TMP_DIR}"

if [ "$APP_MODE" = "dev" ]; then
    exec python "${APP_DIR}/manage.py" runserver "0.0.0.0:${APP_PORT}"
fi

if command -v supercronic >/dev/null 2>&1; then
    # Ensure django-crontab definitions are installed in user crontab first.
    python "${APP_DIR}/manage.py" crontab add >/dev/null 2>&1 || true
    crontab -l 2>/dev/null | sed '/^\s*#/d; /^\s*$/d' > "${CRON_FILE}" || true
    if [ -s "${CRON_FILE}" ]; then
        safe_chmod 600 "${CRON_FILE}"
        : > "${CRON_LOG}"
        supercronic "${CRON_FILE}" > "${CRON_LOG}" 2>&1 &
        CRON_PID=$!
        sleep 1
        if ! kill -0 "${CRON_PID}" 2>/dev/null; then
            echo "supercronic failed to start. Check ${CRON_LOG} for details."
        fi
    else
        echo "No cron entries found. Skipping crond start."
    fi
else
    echo "supercronic not found. Skipping cron."
fi

if [ -n "${WEB_CONCURRENCY:-}" ]; then
    GUNICORN_WORKERS="${WEB_CONCURRENCY}"
else
    cpu_count="$(getconf _NPROCESSORS_ONLN 2>/dev/null || nproc 2>/dev/null || echo 1)"
    if [ "${cpu_count}" -le 2 ]; then
        GUNICORN_WORKERS=2
    else
        GUNICORN_WORKERS=4
    fi
fi

exec gunicorn "${PROJECT_MODULE}.wsgi:application" \
    --bind "0.0.0.0:${APP_PORT}" \
    --workers "${GUNICORN_WORKERS}" \
    --threads "${GUNICORN_THREADS}" \
    --keep-alive "${GUNICORN_KEEPALIVE}" \
    --timeout "${GUNICORN_TIMEOUT}" \
    --worker-tmp-dir /dev/shm
