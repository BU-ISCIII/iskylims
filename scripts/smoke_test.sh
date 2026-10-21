#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"
# Reuse exactly the same engine and Compose frontend selection as the outer
# installer. In particular, Podman prefers podman-compose when it is installed
# instead of delegating `podman compose` to an unrelated Docker Compose plugin.
# shellcheck disable=SC1091
source "$repo_root/deployment/lib/container/common.sh"

install_services=(iskylims)
engine="docker"; mode="production"; compose_file=""; env_file=""
while (($#)); do
    case "$1" in
        --test) mode="test"; shift ;;
        --engine) engine="${2:-}"; shift 2 ;;
        --compose_file) compose_file="${2:-}"; shift 2 ;;
        --env_file) env_file="${2:-}"; shift 2 ;;
        --help) echo "Usage: $0 [--test] [--engine docker|podman] [--compose_file PATH] [--env_file PATH]"; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done
fail() { echo "FAIL: $*" >&2; exit 1; }
compose_file="${compose_file:-$repo_root/docker-compose.$([ "$mode" = test ] && echo test || echo prod).yml}"
env_file="${env_file:-$repo_root/.env.${mode}.file}"
[ -f "$env_file" ] || fail "Compose environment file not found: $env_file. Run container_install.sh first or supply --env_file PATH."
compose_env_file="$env_file"
load_compose_environment_file "$env_file"
set_engine
compose_run() { compose_with_env_exec -f "$compose_file" "$@"; }
compose_run config >/dev/null
    container_id="$(resolve_service_container iskylims)"
    [ -n "$container_id" ] || fail "Service iskylims has no container"
    ensure_service_running iskylims "$container_id" >/dev/null
    engine_exec exec "$container_id" bash -lc 'cd "$INSTALL_PATH" && source virtualenv/bin/activate && python manage.py check && ! python manage.py showmigrations --plan | grep -F '"'"'[ ]'"'"''
    echo "PASS: iskylims Django checks and migrations"
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
