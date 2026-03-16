#!/usr/bin/env bash
set -euo pipefail

APP_DIR="/opt/iskylims"
CRON_DIR="${APP_DIR}/cron"
TMP_DIR="${APP_DIR}/tmp"
CRON_FILE="${CRON_DIR}/iskylims"
CRON_LOG="${TMP_DIR}/supercronic.log"
APP_MODE="${APP_MODE:-prod}"

while [ ! -f "${APP_DIR}/manage.py" ]; do
    sleep 2
done

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
    exec python "${APP_DIR}/manage.py" runserver 0.0.0.0:8001
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

exec gunicorn iskylims.wsgi:application \
    --bind 0.0.0.0:8001 \
    --workers 1 \
    --threads 1 \
    --timeout 120
