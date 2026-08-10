#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"
# Reuse exactly the same engine and Compose frontend selection as the outer
# installer. In particular, Podman prefers podman-compose when it is installed
# instead of delegating `podman compose` to an unrelated Docker Compose plugin.
# shellcheck disable=SC1091
source "$repo_root/deployment/lib/container/common.sh"

install_services=(app)
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
compose_env_file="$env_file"
set_engine
# The generated dotenv file is mode 0600 and contains shell-safe quoted values.
# Source it so direct host checks use the same service ports as Compose.
if [ -n "$env_file" ]; then
    set -a
    # shellcheck disable=SC1090
    source "$env_file"
    set +a
fi
compose_run() { compose_with_env_exec -f "$compose_file" "$@"; }
fail() { echo "FAIL: $*" >&2; exit 1; }
compose_run config >/dev/null
    container_id="$(resolve_service_container app)"
    [ -n "$container_id" ] || fail "Service app has no container"
    ensure_service_running app "$container_id" >/dev/null
    engine_exec exec "$container_id" bash -lc 'cd "$INSTALL_PATH" && source virtualenv/bin/activate && python manage.py check && ! python manage.py showmigrations --plan | grep -F '"'"'[ ]'"'"''
    echo "PASS: app Django checks and migrations"
check_url() {
    local service="$1" url="$2"
    curl --fail --silent --show-error --location --max-time 20 --output /dev/null "$url" \
        || { echo "FAIL: $service health endpoint: $url" >&2; return 1; }
    echo "PASS: $service health endpoint"
}
for service in "${install_services[@]}"; do
    prefix="${service^^}"
    prefix="${prefix//-/_}"
    port_variable="${prefix}_APP_PORT"
    port="${!port_variable:-}"
    [ -n "$port" ] || fail "$port_variable is required in the rendered service settings"
    check_url "$service" "http://127.0.0.1:${port}/health/"
done
echo "iSkyLIMS deployment smoke test passed."
