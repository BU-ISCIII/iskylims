#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "$script_dir/deployment/lib/container/common.sh"
# shellcheck disable=SC1091
source "$script_dir/deployment/lib/container/django.sh"

APP_VERSION="0.1.0"
APPLICATION_NAME="iSkyLIMS"

# ============================================================================
# GENERATED SERVICE/ADD-ON CUSTOMIZATION
# Regenerate these callbacks from the descriptor; keep application-neutral
# lifecycle mechanics below unchanged.
# ============================================================================
install_services=(app)
addon_build_services=()
permission_services=(app apache)
configured_services=(app apache samba)

default_service_install_conf() {
    case "$1" in
        app) [ "$mode" = test ] && echo conf/docker_test_settings.txt || echo conf/docker_production_settings.txt ;;
        apache) [ "$mode" = test ] && echo conf/apache/apache_test_settings.txt || echo conf/apache/apache_production_settings.txt ;;
        samba) [ "$mode" = test ] && echo conf/samba/samba_test_settings.txt || echo conf/samba/samba_production_settings.txt ;;
        *) return 1 ;;
    esac
}
service_build_context_dir() {
    case "$1" in
        app) echo . ;;
        *) return 1 ;;
    esac
}
service_environment_prefix() {
    local prefix
    array_contains "$1" "${install_services[@]}" || return 1
    prefix="${1^^}"
    printf '%s\n' "${prefix//-/_}"
}
service_environment_value() {
    local prefix variable
    prefix="$(service_environment_prefix "$1")" || return 1
    variable="${prefix}_$2"
    if [ -n "${!variable:-}" ]; then
        printf '%s\n' "${!variable}"
    elif [ "$#" -ge 3 ]; then
        printf '%s\n' "$3"
    else
        die "$variable is required in the rendered service settings"
    fi
}
service_repo_path() {
    service_environment_value "$1" REPO_PATH
}
service_install_path() {
    service_environment_value "$1" INSTALL_PATH
}
service_readiness_path() {
    case "$1" in
        app) echo "$(service_install_path "$1")/manage.py" ;;
        *) return 1 ;;
    esac
}
service_image_name() {
    case "$1" in
        app) echo relecov-iskylims:local ;;
        *) return 1 ;;
    esac
}
service_profile() {
    case "$1" in
        app) echo django ;;
        *) return 1 ;;
    esac
}
service_dockerfile() {
    case "$1" in
        app) echo Dockerfile ;;
        *) return 1 ;;
    esac
}
service_container_install_conf() {
    case "$1" in
        app) echo conf/.runtime_install_settings.txt ;;
        *) return 1 ;;
    esac
}
service_uid() {
    service_environment_value "$1" APP_UID
}
service_gid() {
    service_environment_value "$1" APP_GID
}

prepare_compose_environment() {
    local -a settings_sources=(
        "APP|${install_conf_host_by_service[app]}"
        "|${install_conf_host_by_service[apache]}"
        "|${install_conf_host_by_service[samba]}"
    )
    local -a deployment_values=(
        "GIT_REVISION|$git_revision"
        "APP_IMAGE|relecov-iskylims:local"
    )
    compose_env_file="$script_dir/.env.${mode}.file"
    write_compose_environment_file "$compose_env_file" settings_sources deployment_values
}

# Every project uses the generated interpolation file in both modes because it
# combines service-specific settings sources with collision-safe prefixes.
deployment_compose() {
    compose_exec --env-file "$compose_env_file" "$@"
}
current_service_container() {
    resolve_service_container "$1"
}

# Each Django service renders its own protected host settings bind. React
# services and add-ons have no Django settings source.
prepare_application_host_sources() {
    local settings_output
    if [ "$mode" = production ]; then
        settings_output="$(service_environment_value app DJANGO_SETTINGS_PATH)"
        [ -n "$settings_output" ] || { echo "DJANGO_SETTINGS_PATH is required for app" >&2; return 1; }
        mkdir -p "$(dirname "$settings_output")"
        prepare_django_settings_bind_mount ./conf/template_settings.py "$settings_output" "${install_conf_host_by_service[app]}"
    fi
    # conf/apache contains the application-owned Apache sources. Render every
    # deployment value only after the protected settings environment is loaded,
    # then expose the completed files as Compose bind sources.
    local apache_source_dir="$script_dir/conf/apache"
    local apache_output_dir="$script_dir/deployment/apache"
    local apache_conf_name apache_config_service apache_log_path
    apache_config_service=app
    [ -d "$apache_source_dir" ] || {
        echo "Apache source configuration directory not found: $apache_source_dir" >&2
        return 1
    }
    mkdir -p "$apache_output_dir"

    export APACHE_SERVER_NAME="${APACHE_SERVER_NAME:?APACHE_SERVER_NAME is required}"
    export APACHE_UPSTREAM_SERVICE="${APACHE_UPSTREAM_SERVICE:-$apache_config_service}"
    export APACHE_UPSTREAM_PORT="${APACHE_UPSTREAM_PORT:-$(service_environment_value "$apache_config_service" APP_PORT)}"
    # For the default route, INSTALL_PATH means the service selected by
    # ADDONS.apache.CONFIG_SERVICE. Multi-service routes use their explicit
    # API_INSTALL_PATH, WEB_INSTALL_PATH, etc. values instead.
    export INSTALL_PATH="$(service_install_path "$apache_config_service")"
    export APACHE_PROXY_TIMEOUT="${APACHE_PROXY_TIMEOUT:-$(service_environment_value "$apache_config_service" GUNICORN_TIMEOUT 120)}"
    export APACHE_LOG_STEM="${APACHE_LOG_STEM:-$(normalize_apache_server_name "$APACHE_SERVER_NAME")}"

    for apache_conf_name in 00-logs.conf 01-reverse-proxy.conf 02-server-status.conf; do
        [ -f "$apache_source_dir/$apache_conf_name" ] || {
            echo "Apache source configuration not found: $apache_source_dir/$apache_conf_name" >&2
            return 1
        }
        render_environment_config_template \
            "$apache_source_dir/$apache_conf_name" \
            "$apache_output_dir/$apache_conf_name" 0644 || return 1
    done

    # Production bind-mounts Apache logs from the host; tests use a named volume.
    if [ "$mode" = production ]; then
        apache_log_path="${APACHE_LOG_PATH:?APACHE_LOG_PATH is required}"
        mkdir -p "$apache_log_path"
    fi
}

# Keep one independently reviewable host permission specification per
# application and per selected add-on. Empty add-on specs are intentional until
# that add-on declares writable bind sources.
prepare_host_bind_source_permissions() {
    [ "$mode" = production ] || return 0
    local log_path settings_path uid gid
    log_path="$(service_environment_value app HOST_LOG_PATH)"
    settings_path="$(service_environment_value app DJANGO_SETTINGS_PATH)"
    [ -n "$log_path" ] || { echo "HOST_LOG_PATH is required for app" >&2; return 1; }
    [ -n "$settings_path" ] || { echo "DJANGO_SETTINGS_PATH is required for app" >&2; return 1; }
    uid="$(service_uid app)"; gid="$(service_gid app)"
    local -a app_host_bind_permission_spec=(
        "$log_path|$uid:$gid|0775"
        "$(dirname "$settings_path")|-|0755"
        "$settings_path|$uid:$gid|0664"
    )
    apply_host_permission_spec "${app_host_bind_permission_spec[@]}"
    # Generated proxy configuration is read-only in Apache. Its host files need
    # traversal/read permissions, while the production log bind must be writable.
    apache_log_path="${APACHE_LOG_PATH:?APACHE_LOG_PATH is required}"
    local -a apache_host_bind_permission_spec=(
        "$script_dir/deployment/apache|-|0755"
        "$script_dir/deployment/apache/00-logs.conf|-|0644"
        "$script_dir/deployment/apache/01-reverse-proxy.conf|-|0644"
        "$script_dir/deployment/apache/02-server-status.conf|-|0644"
        # registry.access.redhat.com/ubi9/httpd-24 runs as UID 1001 with GID 0.
        # The shared helper applies these IDs directly for Docker and through
        # podman unshare when the bind source belongs to a rootless userns.
        "$apache_log_path|1001:0|0775"
    )
    apply_host_permission_spec "${apache_host_bind_permission_spec[@]}"
}

# Keep a separate running-mount specification in every service/add-on case.
prepare_running_container_mount_permissions() {
    local service_name="$1" container_id="$2"
    local install_path uid gid
    case "$service_name" in
        app)
            install_path="$(service_install_path "$service_name")"
            uid="$(service_uid "$service_name")"; gid="$(service_gid "$service_name")"
            local -a app_running_mount_permission_spec=(
                "$install_path/logs|$uid:$gid|u+rwX,g+rwX"
                "$install_path/documents|$uid:$gid|u+rwX,g+rwX"
                "$install_path/static|$uid:$gid|u+rwX,g+rwX,o+rX"
            )
            apply_container_directory_permission_spec "$container_id" "${app_running_mount_permission_spec[@]}"
            prepare_django_container_settings_permissions "$container_id" "$install_path/iskylims/settings.py" "$uid" "$gid"
            ;;
        apache)
            # Apache currently needs no ownership repair inside its running
            # container. Keep an explicit add-on policy ready for future mounts.
            local -a apache_running_mount_permission_spec=()
            apply_container_directory_permission_spec "$container_id" "${apache_running_mount_permission_spec[@]}"
            ;;
        *) return 0 ;;
    esac
}

bootstrap_service() {
    local service_name="$1" container_id="$2" deployment_action="$3"
    local repo_path runtime_conf uid gid status
    local -a args
    case "$service_name" in
        app)
            repo_path="$(service_repo_path "$service_name")"
            # Fixed temporary in-container path; this is not operator configuration.
            runtime_conf=conf/.runtime_install_settings.txt
            [[ "$runtime_conf" == /* ]] || runtime_conf="$repo_path/$runtime_conf"
            uid="$(service_uid "$service_name")"; gid="$(service_gid "$service_name")"
            stage_container_runtime_config "$container_id" "${install_conf_host_by_service[$service_name]}" "$runtime_conf" "$uid" "$gid"
            args=(--bootstrap "$deployment_action" --git_revision "$git_revision" --conf "$runtime_conf" --skip_apache_restart)
            for hook in "${migration_script_before[@]}"; do args+=(--script_before "$hook"); done
            for hook in "${migration_script_after[@]}"; do args+=(--script_after "$hook"); done
            status=0; engine_exec exec "$container_id" bash "$repo_path/install.sh" "${args[@]}" || status=$?
            [ "$mode" = test ] || remove_container_runtime_config "$container_id" "$runtime_conf" || true
            return "$status"
            ;;
        *) return 0 ;;
    esac
}

# Applications with disposable fixtures or demo files customize this callback
# in their generated wrapper and set application_supports_test_data=true. Keep
# application-specific fixture names, users/groups, downloads, and data-service
# layout here so the complete test installation remains readable in one file.
application_supports_test_data=true
load_test_deployment_data() {
    local app_container app_install_path samba_container archive downloaded_archive
    local admin_groups_code

    app_container="$(current_service_container app)" \
        || die "Unable to resolve the app container for test-data loading"
    app_install_path="$(service_install_path app)"

    if [ "$skip_test_data" = false ]; then
        echo "Loading iSkyLIMS test fixtures"
        engine_exec exec -w "$app_install_path" "$app_container" \
            "$app_install_path/virtualenv/bin/python" manage.py \
            loaddata test/test_data.json
        admin_groups_code=$(cat <<'PY'
from django.contrib.auth.models import Group, User

admin = User.objects.get(username="admin")
admin.groups.add(
    Group.objects.get(name="WetlabManager"),
    Group.objects.get(name="ServiceManager"),
)
print("admin groups:", list(admin.groups.values_list("name", flat=True)))
PY
)
        engine_exec exec -w "$app_install_path" "$app_container" \
            "$app_install_path/virtualenv/bin/python" manage.py \
            shell -c "$admin_groups_code"
    else
        echo "Skipping iSkyLIMS test fixtures as requested"
    fi

    if [ "$skip_demo_data" = true ]; then
        echo "Skipping Samba demo-data load as requested"
        return 0
    fi

    samba_container="$(current_service_container samba)" \
        || die "The iSkyLIMS test-data workflow requires the Samba service"
    archive="$demo_data"
    downloaded_archive=""
    if [ -z "$archive" ]; then
        command -v wget >/dev/null 2>&1 \
            || die "wget is required to download iSkyLIMS demo data"
        downloaded_archive="$(mktemp /tmp/iskylims_demo_data.XXXXXX.tar.gz)"
        archive="$downloaded_archive"
        wget -O "$archive" \
            https://zenodo.org/record/8091169/files/iskylims_demo_data.tar.gz
    fi
    [ -f "$archive" ] || die "Demo-data archive not found: $archive"

    echo "Copying and extracting iSkyLIMS demo data in Samba"
    engine_exec cp "$archive" "$samba_container:/mnt/iskylims_demo_data.tar.gz"
    engine_exec exec "$samba_container" \
        tar -xf /mnt/iskylims_demo_data.tar.gz -C /mnt
    engine_exec exec "$samba_container" sh -lc '
        for root in /mnt/test_ngs_data /mnt/Runs; do
            if [ -d "$root" ]; then
                find "$root" -type d -exec chmod o+rx {} +
                find "$root" -type f -exec chmod o+r {} +
            fi
        done
    '
    engine_exec exec "$samba_container" rm /mnt/iskylims_demo_data.tar.gz
    [ -z "$downloaded_archive" ] || rm -f "$downloaded_archive"
}

action="install"; mode="production"; engine="docker"; git_revision="current"
install_conf=""; compose_file=""; compose_env_file=""
install_conf_map_entries=(); migration_script_before=(); migration_script_after=()
demo_data=""; skip_demo_data=""; skip_test_data=""

usage() {
    cat <<'EOF'
Install, upgrade, or repair the application deployment.

Options:
  --action install|upgrade|fix-permissions
  --test
  --engine docker|podman
  --git_revision <branch|tag|commit|current>
  --install_conf <path>              First application service only.
  --install_conf_map <component,path>  Repeat for application and add-on overrides.
  --compose_file <path>
  --script_before <name[,args]>
  --script_after <name[,args]>
  --script <name[,args]>
  --demo_data <path>                 Import application demo data on install.
  --skip_demo_data
  --skip_test_data
  --help
  --version
EOF
}
die() { echo "ERROR: $*" >&2; exit 1; }

# 1. Parse the canonical outer-installer interface.
while (($#)); do
    case "$1" in
        --action) action="${2:-}"; shift 2 ;;
        --test) mode="test"; shift ;;
        --engine) engine="${2:-}"; shift 2 ;;
        --git_revision) git_revision="${2:-}"; shift 2 ;;
        --install_conf) install_conf="${2:-}"; shift 2 ;;
        --install_conf_map) install_conf_map_entries+=("${2:-}"); shift 2 ;;
        --compose_file) compose_file="${2:-}"; shift 2 ;;
        --script_before) migration_script_before+=("${2:-}"); shift 2 ;;
        --script_after|--script) migration_script_after+=("${2:-}"); shift 2 ;;
        --demo_data) demo_data="${2:-}"; shift 2 ;;
        --skip_demo_data) skip_demo_data=true; shift ;;
        --skip_test_data) skip_test_data=true; shift ;;
        --help) usage; exit 0 ;;
        --version) echo "$APP_VERSION"; exit 0 ;;
        *) die "Unknown option: $1" ;;
    esac
done

# 2. Validate arguments before modifying deployment state.
[[ "$action" =~ ^(install|upgrade|fix-permissions)$ ]] || die "Invalid action: $action"
[[ "$engine" =~ ^(docker|podman)$ ]] || die "Invalid engine: $engine"
if [ -n "$demo_data" ] && [ "$application_supports_test_data" != true ]; then
    die "--demo_data is not implemented for $APPLICATION_NAME"
fi
if [ -n "$demo_data" ]; then
    [ -f "$demo_data" ] || die "Demo-data file not found: $demo_data"
    demo_data="$(cd "$(dirname "$demo_data")" && pwd)/$(basename "$demo_data")"
fi
if [ "$mode" = test ] && [ "$action" = install ] \
    && [ "$application_supports_test_data" = true ]; then
    skip_demo_data="${skip_demo_data:-false}"
    skip_test_data="${skip_test_data:-false}"
elif [ "$action" = install ] && [ -n "$demo_data" ] \
    && [ "$application_supports_test_data" = true ]; then
    skip_demo_data="${skip_demo_data:-false}"
    skip_test_data=true
else
    skip_demo_data=true
    skip_test_data=true
fi

# 3. Resolve one protected configuration source per configured component.
cd "$script_dir"
declare -A install_conf_host_by_service=()
for service_name in "${configured_services[@]}"; do
    install_conf_host_by_service["$service_name"]="$(default_service_install_conf "$service_name")"
done
if [ -n "$install_conf" ]; then install_conf_host_by_service["${install_services[0]}"]="$install_conf"; fi
for mapping in "${install_conf_map_entries[@]}"; do
    [[ "$mapping" == *,* ]] || die "Invalid --install_conf_map: $mapping"
    service_name="${mapping%%,*}"; path="${mapping#*,}"
    array_contains "$service_name" "${configured_services[@]}" || die "Unknown mapped component: $service_name"
    install_conf_host_by_service["$service_name"]="$path"
done
for service_name in "${configured_services[@]}"; do
    path="${install_conf_host_by_service[$service_name]}"
    [[ "$path" = /* ]] || path="$script_dir/$path"
    [ -f "$path" ] || die "Configuration for $service_name not found: $path"
    if [ "$mode" = production ] \
        && grep -Eq '^[A-Z0-9_]+=.*CHANGE_ME' "$path"; then
        die "Production configuration for $service_name contains CHANGE_ME: $path"
    fi
    install_conf_host_by_service["$service_name"]="$(cd "$(dirname "$path")" && pwd)/$(basename "$path")"
done

# 4. Select the engine, prepare host sources and validate the final Compose model.
set_engine "$engine"
compose_file="${compose_file:-docker-compose.$([ "$mode" = test ] && echo test || echo prod).yml}"
require_compose_file "$compose_file"
prepare_compose_environment
load_compose_environment_file "$compose_env_file"
prepare_application_host_sources
prepare_host_bind_source_permissions
deployment_compose -f "$compose_file" config \
    || die "Compose configuration validation failed: $compose_file"

# 5. Dispatch permission-only repair without building or bootstrapping.
if [ "$action" = fix-permissions ]; then
    for service_name in "${permission_services[@]}"; do
        container_id="$(current_service_container "$service_name" 2>/dev/null || true)"
        [ -z "$container_id" ] || prepare_running_container_mount_permissions "$service_name" "$container_id"
    done
    echo "Permissions repaired without build or bootstrap."
    exit 0
fi

# 6. Build application services in declared order. Production builds use the
# engine directly: Django receives its settings as
# an ephemeral build secret, while React receives only its public VITE value.
# This avoids requiring Compose implementations to support build.secrets.
for service_name in "${install_services[@]}"; do
    if [ "$mode" = test ]; then
        deployment_compose -f "$compose_file" build --no-cache "$service_name"
        continue
    fi
    context="$(service_build_context_dir "$service_name")"
    dockerfile="$(service_dockerfile "$service_name")"
    profile="$(service_profile "$service_name")"
    if [ "$profile" = django ]; then
        engine_build --no-cache --file "$context/$dockerfile" \
            --secret "id=install_conf,src=${install_conf_host_by_service[$service_name]}" \
            --build-arg GIT_REVISION="$git_revision" \
            --build-arg INSTALL_CONF="$(service_container_install_conf "$service_name")" \
            --build-arg USE_INSTALL_CONF_SECRET=true \
            --build-arg RENDER_DJANGO_SETTINGS=false \
            --build-arg APP_REPO_PATH="$(service_repo_path "$service_name")" \
            --build-arg APP_INSTALL_PATH="$(service_install_path "$service_name")" \
            --build-arg APP_PORT="$(service_environment_value "$service_name" APP_PORT)" \
            --build-arg APP_UID="$(service_uid "$service_name")" \
            --build-arg APP_GID="$(service_gid "$service_name")" \
            --tag "$(service_image_name "$service_name")" "$context"
    else
        vite_api_url="$(service_environment_value "$service_name" VITE_API_BASE_URL)"
        engine_build --no-cache --file "$context/$dockerfile" \
            --build-arg GIT_REVISION="$git_revision" \
            --build-arg VITE_API_BASE_URL="$vite_api_url" \
            --tag "$(service_image_name "$service_name")" "$context"
    fi
done
# Build add-on images through Compose so their declared build arguments and
# add-on-owned Dockerfiles remain the single source of truth.
for service_name in "${addon_build_services[@]}"; do
    deployment_compose -f "$compose_file" build --no-cache "$service_name"
done
# 7. Recreate and start the complete topology from one Compose invocation so
# freshly built images and the current configuration are deployed consistently.
# Named volumes and bind-mounted persistent data are preserved.
deployment_compose -f "$compose_file" up -d --force-recreate

# 8. Wait for every application service readiness contract.
for service_name in "${install_services[@]}"; do
    container_id="$(current_service_container "$service_name")"
    [ -n "$container_id" ] || die "Unable to resolve $service_name container"
    ensure_service_running "$service_name" "$container_id" >/dev/null
    readiness_path="$(service_readiness_path "$service_name")"
    deadline=$((SECONDS + 120))
    until engine_exec exec "$container_id" test -f "$readiness_path"; do
        ((SECONDS < deadline)) || { engine_exec logs --tail 200 "$container_id"; die "$service_name readiness timeout"; }
        sleep 2
    done
done

# 9. Repair running mounts for applications and selected add-ons.
for service_name in "${permission_services[@]}"; do
    container_id="$(current_service_container "$service_name")"
    prepare_running_container_mount_permissions "$service_name" "$container_id"
done

# 10. Bootstrap only application profiles that require runtime bootstrap.
for service_name in "${install_services[@]}"; do
    container_id="$(current_service_container "$service_name")"
    bootstrap_service "$service_name" "$container_id" "$action" || die "$service_name bootstrap failed"
done

# 11. Load application-owned data for a fresh test install, or for an explicit
# production --demo_data request. Production never imports data implicitly.
if [ "$action" = install ] \
    && [ "$application_supports_test_data" = true ] \
    && { [ "$mode" = test ] || [ -n "$demo_data" ]; }; then
    load_test_deployment_data
fi

# 12. Execute the common smoke dispatcher with generated profile checks.
smoke_args=(--engine "$engine" --compose_file "$compose_file" --env_file "$compose_env_file")
[ "$mode" = test ] && smoke_args+=(--test)
bash "$script_dir/scripts/smoke_test.sh" "${smoke_args[@]}"
echo "$action completed successfully for $APPLICATION_NAME."
