#!/usr/bin/env bash
set -euo pipefail

engine="docker"; mode="production"; compose_file=""; env_file=""
while (($#)); do
    case "$1" in
        --test) mode="test"; shift ;;
        --engine) engine="${2:-}"; shift 2 ;;
        --compose_file) compose_file="${2:-}"; shift 2 ;;
        --env_file) env_file="${2:-}"; shift 2 ;;
        --help) echo "Usage: $0 [--test] [--engine docker|podman]"; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done
compose_file="${compose_file:-docker-compose.$([ "$mode" = test ] && echo test || echo prod).yml}"
if [ "$engine" = docker ]; then compose=(docker compose); else compose=(podman compose); fi
args=(-f "$compose_file"); [ -z "$env_file" ] || args=(--env-file "$env_file" "${args[@]}")
# The generated dotenv file is mode 0600 and contains shell-safe quoted values.
# Source it so direct host checks use the same service ports as Compose.
if [ -n "$env_file" ]; then
    set -a
    # shellcheck disable=SC1090
    source "$env_file"
    set +a
fi
compose_run() { "${compose[@]}" "${args[@]}" "$@"; }
fail() { echo "FAIL: $*" >&2; exit 1; }
compose_run config --quiet
    container_id="$(compose_run ps -q app)"
    [ -n "$container_id" ] || fail "Service app has no container"
    compose_run exec -T app bash -lc 'cd "$INSTALL_PATH" && source virtualenv/bin/activate && python manage.py check && ! python manage.py showmigrations --plan | grep -F '"'"'[ ]'"'"''
    echo "PASS: app Django checks and migrations"
check_url() {
    local service="$1" url="$2"
    curl --fail --silent --show-error --location --max-time 20 --output /dev/null "$url" \
        || { echo "FAIL: $service health endpoint: $url" >&2; return 1; }
    echo "PASS: $service health endpoint"
}
    check_url app "http://127.0.0.1:${APP_APP_PORT:-8000}/health/"
echo "iSkyLIMS deployment smoke test passed."
