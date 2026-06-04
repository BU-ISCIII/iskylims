#!/usr/bin/env bash
set -euo pipefail

APP_DIR="${INSTALL_PATH:-/opt/iskylims}"
CRON_DIR="${APP_DIR}/cron"
TMP_DIR="${APP_DIR}/tmp"
CRON_FILE="${CRON_DIR}/iskylims"
CRON_LOG="${TMP_DIR}/supercronic.log"
CRON_DISABLED_FILE="${CRON_DIR}/disabled"
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

if [ -f "${CRON_DISABLED_FILE}" ]; then
    echo "Cron is disabled by ${CRON_DISABLED_FILE}. Skipping supercronic start."
elif command -v supercronic >/dev/null 2>&1; then
    # Build supercronic's crontab directly from Django settings. Avoid the
    # system crontab command because it depends on PAM behavior that varies
    # across rootless container hosts.
    python - <<'PY' > "${CRON_FILE}"
import os
import shlex

os.environ.setdefault(
    "DJANGO_SETTINGS_MODULE",
    os.environ.get("DJANGO_SETTINGS_MODULE", "iskylims.settings"),
)

import django
django.setup()

from django.conf import settings

app_dir = os.environ.get("INSTALL_PATH", "/opt/iskylims")
python_bin = os.path.join(app_dir, "virtualenv", "bin", "python")
settings_module = os.environ.get("DJANGO_SETTINGS_MODULE", "iskylims.settings")
command_suffix = getattr(settings, "CRONTAB_COMMAND_SUFFIX", "")

for job in getattr(settings, "CRONJOBS", []):
    if len(job) < 2:
        continue

    schedule = job[0]
    dotted_path = job[1]
    job_suffix = job[2] if len(job) > 2 else ""
    module_name, function_name = dotted_path.rsplit(".", 1)
    python_code = (
        "import os; "
        f"os.environ.setdefault('DJANGO_SETTINGS_MODULE', {settings_module!r}); "
        "import django; django.setup(); "
        f"from {module_name} import {function_name} as cron_job; "
        "cron_job()"
    )
    command = (
        f"cd {shlex.quote(app_dir)} && "
        f"DJANGO_SETTINGS_MODULE={shlex.quote(settings_module)} "
        f"{shlex.quote(python_bin)} -c {shlex.quote(python_code)}"
    )
    suffixes = " ".join(s for s in (job_suffix, command_suffix) if s)
    if suffixes:
        command = f"{command} {suffixes}"

    print(f"{schedule} {command}")
PY
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
