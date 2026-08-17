#!/usr/bin/env bash

# BU-ISCIII shared host-side Django deployment helpers.
#
# These functions belong to container_install.sh, not the in-container
# install.sh: Compose requires the settings bind-mount source to exist on the
# host before the application container can start. The inner installer cannot
# create or repair that host path because it runs after mounts are established.

BU_ISCIII_DJANGO_CONTAINER_LIB_VERSION="0.2.0"

# Generate a cryptographically random Django-compatible SECRET_KEY without
# requiring Django itself to be installed on the deployment host. It belongs in
# this profile because its alphabet and purpose are Django-specific.
generate_django_secret_key() {
    if command -v python3 >/dev/null 2>&1; then
        python3 -c "import secrets; print(''.join(secrets.choice('abcdefghijklmnopqrstuvwxyz0123456789!@#$%^&*(-_=+)') for _ in range(50)))"
    else
        LC_ALL=C tr -dc 'A-Za-z0-9!@#$%^&*(-_=+)' < /dev/urandom | head -c 50
        printf "\n"
    fi
}

# Render a Django settings template from the standard BU-ISCIII installation
# keys while preserving an existing non-placeholder SECRET_KEY.
# Arguments: template path, destination settings path, install-conf path.
render_django_settings_file() {
    local template_path="$1"
    local settings_path="$2"
    local install_conf_path="$3"
    local secret_line=""
    local template_secret_line=""
    local django_debug=""
    local token variable value python_literal
    local -a application_replacements=()

    if [ -f "$settings_path" ]; then
        secret_line="$(grep -E "^SECRET_KEY[[:space:]]*=" "$settings_path" | tail -n 1 || true)"
    fi
    if [ -z "$secret_line" ] || [[ "$secret_line" =~ SECRET_KEY[[:space:]]*=[[:space:]]*SECRET ]]; then
        secret_line="SECRET_KEY = '$(generate_django_secret_key)'"
    fi

    template_secret_line="$(grep -E "^SECRET_KEY[[:space:]]*=" "$template_path" | head -n 1)"
    if [ -z "$template_secret_line" ]; then
        echo "Django template '$template_path' has no SECRET_KEY assignment." >&2
        return 1
    fi
    django_debug="$(config_value_or_default DJANGO_DEBUG "$install_conf_path" false)"
    case "${django_debug,,}" in
        true|1|yes|on) django_debug=True ;;
        false|0|no|off) django_debug=False ;;
        *)
            echo "DJANGO_DEBUG must be a boolean value." >&2
            return 1
            ;;
    esac
    while IFS= read -r token; do
        [ -n "$token" ] || continue
        variable="${token#settingsconf_}"
        value="$(read_install_conf_value "$variable" "$install_conf_path")" || return 1
        python_literal="$(python3 -c \
            'import json, sys; print(json.dumps(sys.argv[1]))' "$value")" || return 1
        application_replacements+=("$token" "$python_literal")
    done < <(grep -oE 'settingsconf_[A-Z][A-Z0-9_]*' "$template_path" | sort -u || true)
    render_config_template "$template_path" "$settings_path" 0664 \
        "$template_secret_line" "$secret_line" \
        djangouser "$(read_install_conf_first "$install_conf_path" DB_USER)" \
        djangopass "$(read_install_conf_value DB_PASSWORD "$install_conf_path")" \
        djangohost "$(read_install_conf_value DB_HOST "$install_conf_path")" \
        djangoport "$(read_install_conf_value DB_PORT "$install_conf_path")" \
        djangodbname "$(read_install_conf_value DB_NAME "$install_conf_path")" \
        emailhostserver "$(read_install_conf_value EMAIL_HOST "$install_conf_path")" \
        emailport "$(read_install_conf_value EMAIL_PORT "$install_conf_path")" \
        emailhostuser "$(read_install_conf_value EMAIL_HOST_USER "$install_conf_path")" \
        emailhostpassword "$(read_install_conf_value EMAIL_HOST_PASSWORD "$install_conf_path")" \
        emailhosttls "$(read_install_conf_value EMAIL_USE_TLS "$install_conf_path")" \
        djangodebug "$django_debug" \
        djangoallowedhosts "$(read_install_conf_value DJANGO_ALLOWED_HOSTS "$install_conf_path")" \
        djangocsrftrustedorigins "$(read_install_conf_value DJANGO_CSRF_TRUSTED_ORIGINS "$install_conf_path")" \
        dbconnmaxage "$(config_value_or_default DB_CONN_MAX_AGE "$install_conf_path" 0)" \
        "${application_replacements[@]}"
}

# Ensure the production settings bind source exists and reflects the selected
# database configuration before Compose starts the application container.
# Arguments: template path, destination settings path, install-conf path.
prepare_django_settings_bind_mount() {
    local template_path="$1"
    local settings_path="$2"
    local install_conf_path="$3"

    [ "${mode:-production}" = "production" ] || return 0
    if [ -d "$settings_path" ]; then
        echo "DJANGO_SETTINGS_PATH must resolve to a file path, but '$settings_path' is a directory." >&2
        return 1
    fi

    mkdir -p "$(dirname "$settings_path")"
    # Always rerender so application-template and deployment-setting changes
    # reach upgrades. render_django_settings_file preserves the existing
    # non-placeholder SECRET_KEY, and render_config_template installs atomically.
    render_django_settings_file "$template_path" "$settings_path" "$install_conf_path"
    chmod_with_engine_fallback 0664 "$settings_path"
}

# Apply ownership and mode to a Django settings file visible inside a running
# container. The caller supplies the exact application path and identity;
# writable-directory policy remains application-owned.
# Arguments: container ID/name, settings path, application UID, application GID.
prepare_django_container_settings_permissions() {
    local container_id="$1"
    local settings_path="$2"
    local application_uid="$3"
    local application_gid="$4"

    [ -n "$settings_path" ] || return 0
    engine_exec exec --user 0 "$container_id" sh -c '
        settings_path="$1"
        owner="$2"
        if [ -f "$settings_path" ]; then
            chown "$owner" "$settings_path"
            chmod 0664 "$settings_path"
        fi
    ' _ "$settings_path" "$application_uid:$application_gid"
}
