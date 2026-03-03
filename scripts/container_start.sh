#!/usr/bin/env bash
set -euo pipefail

APP_DIR="/opt/iskylims"
CRON_DIR="${APP_DIR}/cron"
TMP_DIR="${APP_DIR}/tmp"
CRON_FILE="${CRON_DIR}/iskylims"
CRON_LOG="${TMP_DIR}/crond.log"
APP_MODE="${APP_MODE:-prod}"

while [ ! -f "${APP_DIR}/manage.py" ]; do
    sleep 2
done

source "${APP_DIR}/virtualenv/bin/activate"

mkdir -p "${CRON_DIR}" "${TMP_DIR}"
chmod 700 "${CRON_DIR}" "${TMP_DIR}"

if [ "$APP_MODE" = "dev" ]; then
    exec python "${APP_DIR}/manage.py" runserver 0.0.0.0:8001
fi

if command -v crond >/dev/null 2>&1; then
    # django-crontab needs an explicit add before jobs are visible via `crontab show`
    python "${APP_DIR}/manage.py" crontab add >/dev/null 2>&1 || true
    python "${APP_DIR}/manage.py" crontab show 2>/dev/null \
        | sed '/^no crontab for /d; /^\s*$/d' > "${CRON_FILE}" || true
    if [ -s "${CRON_FILE}" ]; then
        chmod 600 "${CRON_FILE}"
        : > "${CRON_LOG}"
        crond -n -m off -c "${CRON_DIR}" > "${CRON_LOG}" 2>&1 &
        CROND_PID=$!
        sleep 1
        if ! kill -0 "${CROND_PID}" 2>/dev/null; then
            if grep -qi "permission denied" "${CRON_LOG}" 2>/dev/null; then
                echo "crond could not start under non-root user. Skipping cron daemon startup."
            else
                echo "crond failed to start. Check ${CRON_LOG} for details."
            fi
        fi
    else
        echo "No cron entries found. Skipping crond start."
    fi
else
    echo "crond not found. Skipping cron."
fi

exec gunicorn iskylims.wsgi:application \
    --bind 0.0.0.0:8001 \
    --workers 3 \
    --threads 2 \
    --timeout 120
