#!/usr/bin/env bash

# BU-ISCIII shared container installer library.
#
# This file is centrally managed. Applications must not edit vendored copies.
# Application container_install.sh scripts provide the global variables used by
# these helpers (engine, mode, compose_env_file, ENGINE_CMD and COMPOSE_CMD).

BU_ISCIII_CONTAINER_LIB_VERSION="0.1.0"

# Select and validate the requested container engine and its Compose frontend.
# This keeps Docker/Podman detection identical in every outer installer.
set_engine() {
    if [ "$engine" = "docker" ]; then
        if ! command -v docker >/dev/null 2>&1; then
            echo "docker not found. Install docker or use --engine podman."
            exit 1
        fi
        ENGINE_CMD=("docker")
        COMPOSE_CMD=("docker" "compose")
    else
        if ! command -v podman >/dev/null 2>&1; then
            echo "podman not found. Install podman or use --engine docker."
            exit 1
        fi
        ENGINE_CMD=("podman")
        if command -v podman-compose >/dev/null 2>&1; then
            COMPOSE_CMD=("podman-compose")
        elif podman compose version >/dev/null 2>&1; then
            COMPOSE_CMD=("podman" "compose")
        else
            echo "podman compose not available. Install podman-compose or use --engine docker."
            exit 1
        fi
    fi
}

# Run a container-engine command using the engine selected by set_engine.
engine_exec() {
    "${ENGINE_CMD[@]}" "$@"
}

# Build an image with the selected engine. Docker explicitly enables BuildKit
# because Dockerfile secret mounts are unavailable in the legacy builder.
engine_build() {
    if [ "${engine:-docker}" = "docker" ]; then
        DOCKER_BUILDKIT=1 engine_exec build "$@"
    else
        engine_exec build "$@"
    fi
}

# Run a Compose command using the frontend selected by set_engine.
compose_exec() {
    "${COMPOSE_CMD[@]}" "$@"
}

# Run Compose with the generated production environment file when available.
# Test Compose files intentionally continue to use their normal environment.
compose_with_env_exec() {
    if [ "${mode:-production}" = "production" ] \
        && [ -n "${compose_env_file:-}" ] \
        && [ -f "$compose_env_file" ]; then
        compose_exec --env-file "$compose_env_file" "$@"
    else
        compose_exec "$@"
    fi
}

# Return the image ID associated with one service in an explicit Compose file.
# Both arguments are required so the helper has no hidden application globals.
# Arguments: compose file path, service name.
compose_service_image_id() {
    local compose_path="$1"
    local service_name="$2"
    compose_with_env_exec -f "$compose_path" images -q "$service_name" 2>/dev/null \
        | tail -n 1
}

# Fail early when the selected Compose file does not exist.
# Arguments: Compose file path.
require_compose_file() {
    local compose_path="$1"
    if [ ! -f "$compose_path" ]; then
        echo "Compose file '$compose_path' not found" >&2
        return 1
    fi
}

# Validate the fully interpolated Compose model after application environment
# preparation. Arguments: Compose file path.
validate_compose_configuration() {
    local compose_path="$1"
    require_compose_file "$compose_path" || return 1
    if ! compose_with_env_exec -f "$compose_path" config --quiet; then
        echo "Compose configuration validation failed: $compose_path" >&2
        return 1
    fi
}

# Return a repository's full or short HEAD without printing diagnostics.
# Arguments: repository path, optional format (`full` or `short`).
repository_revision() {
    local repository_path="$1"
    local format="${2:-full}"
    command -v git >/dev/null 2>&1 || return 1
    git -C "$repository_path" rev-parse --is-inside-work-tree >/dev/null 2>&1 \
        || return 1
    if [ "$format" = "short" ]; then
        git -C "$repository_path" rev-parse --short HEAD
    else
        git -C "$repository_path" rev-parse HEAD
    fi
}

# Print the actual local checkout and the revision requested by the operator.
# Arguments: display label, repository path, expected branch/tag/commit/current.
print_repository_diagnostics() {
    local label="$1"
    local repository_path="$2"
    local expected_revision="$3"
    local actual_hash=""

    echo "$label"
    echo "  requested revision: $expected_revision"
    actual_hash="$(repository_revision "$repository_path" full)" || true
    if [ -n "$actual_hash" ]; then
        echo "  local HEAD: $(git -C "$repository_path" log -1 --oneline)"
        echo "  local HEAD hash: $actual_hash"
    else
        echo "  local git metadata unavailable"
    fi
}

# Print an image ID before a build. The caller resolves the ID from an explicit
# image name or Compose service and passes it in, keeping lookup policy separate.
# Arguments: display label, image ID/name (empty means no existing image).
print_image_before_diagnostics() {
    local label="$1"
    local image_reference="$2"
    echo "$label"
    if [ -n "$image_reference" ]; then
        echo "  image before build: $image_reference"
    else
        echo "  image before build: none"
    fi
}

# Compare explicit pre-build and post-build image IDs/names.
# Arguments: display label, previous image reference, current image reference.
print_image_after_diagnostics() {
    local label="$1"
    local previous_reference="$2"
    local current_reference="$3"
    echo "$label"
    if [ -z "$current_reference" ]; then
        echo "  image after build: not found"
    elif [ -n "$previous_reference" ] && [ "$previous_reference" = "$current_reference" ]; then
        echo "  image after build: $current_reference"
        echo "  image id check: unchanged"
    elif [ -n "$previous_reference" ]; then
        echo "  image after build: $current_reference"
        echo "  image id check: changed"
    else
        echo "  image after build: $current_reference"
        echo "  image id check: created"
    fi
}

# Print the standard repository and existing-image report before a build.
# Arguments: label, repository path, requested revision, existing image ID/name.
print_prebuild_diagnostics() {
    local label="$1"
    local repository_path="$2"
    local expected_revision="$3"
    local image_reference="$4"

    print_repository_diagnostics "$label source:" "$repository_path" "$expected_revision"
    print_image_before_diagnostics "$label image before build:" "$image_reference"
}

# Compare a container's source checkout with an explicit expected revision.
# Arguments: label, container ID/name, repository path inside the container,
# expected full hash, expected short hash.
print_container_repository_diagnostics() {
    local label="$1"
    local container_id="$2"
    local repository_path="$3"
    local expected_hash="$4"
    local expected_short="$5"
    local container_hash=""
    local container_short=""

    echo "$label"
    container_hash="$(engine_exec exec "$container_id" sh -lc "
        if [ -d '$repository_path/.git' ]; then
            cd '$repository_path' && git rev-parse HEAD
        fi
    " 2>/dev/null | tail -n 1)"
    container_short="$(engine_exec exec "$container_id" sh -lc "
        if [ -d '$repository_path/.git' ]; then
            cd '$repository_path' && git rev-parse --short HEAD
        fi
    " 2>/dev/null | tail -n 1)"
    if [ -n "$container_hash" ]; then
        echo "  container source HEAD hash: $container_hash"
    fi
    if [ -n "$expected_hash" ] && [ -n "$container_hash" ]; then
        if [ "$expected_hash" = "$container_hash" ]; then
            echo "  HEAD check: OK local=$expected_short container=$container_short"
        else
            echo "  HEAD check: MISMATCH local=$expected_short container=$container_short"
        fi
    fi
    engine_exec exec "$container_id" sh -lc "
        echo '  $repository_path HEAD:'
        if [ -d '$repository_path/.git' ]; then
            cd '$repository_path' && git log -1 --oneline
        else
            echo 'not a git checkout'
        fi
    " || true
}

# Copy a host file, entering Podman's user namespace when mapped ownership
# prevents an ordinary host-side copy.
copy_with_podman_fallback() {
    local src="$1"
    local dst="$2"
    local tmp_dst=""

    # Replace the destination through its parent directory instead of opening
    # an existing bind-mounted file in place. Containers may have changed that
    # file's ownership to an unmapped/root identity, while the operator still
    # owns the parent directory and is therefore allowed to replace it. The
    # same-directory rename also prevents readers from seeing a partial file.
    if tmp_dst="$(mktemp "${dst}.tmp.XXXXXX" 2>/dev/null)" \
        && cp "$src" "$tmp_dst" 2>/dev/null \
        && mv -f "$tmp_dst" "$dst" 2>/dev/null; then
        return 0
    fi
    [ -z "$tmp_dst" ] || rm -f "$tmp_dst" 2>/dev/null || true
    if [ "$engine" = "podman" ] && podman unshare cp "$src" "$dst"; then
        return 0
    fi
    echo "Failed to copy '$src' to '$dst'" >&2
    return 1
}

# Run one ownership/mode operation through rootful Docker. Every bind source is
# resolved first and the host root is always rejected.
docker_host_path_operation() {
    local operation="$1"
    local value="$2"
    shift 2
    local path resolved_path

    [[ "$operation" =~ ^(chown|chmod)$ ]] || return 2
    for path in "$@"; do
        resolved_path="$(readlink -f -- "$path" 2>/dev/null || true)"
        if [ -z "$resolved_path" ] || [ "$resolved_path" = / ]; then
            echo "Refusing Docker $operation fallback for unsafe path: $path" >&2
            return 1
        fi
        engine_exec run --rm --user 0 \
            --volume "$resolved_path:/target:z" \
            --entrypoint "/usr/bin/$operation" \
            registry.access.redhat.com/ubi9/ubi-minimal:latest \
            -R "$value" /target || return 1
    done
}

# Change host-path modes, entering the selected engine's ownership context when
# an ordinary host operation is not permitted.
chmod_with_engine_fallback() {
    local mode_value="$1"
    shift

    if chmod "$mode_value" "$@" 2>/dev/null; then
        return 0
    fi
    if [ "$engine" = "podman" ] && podman unshare chmod "$mode_value" "$@"; then
        return 0
    fi
    if [ "$engine" = "docker" ]; then
        [[ "$mode_value" =~ ^[0-7]{3,4}$ ]] || {
            echo "Docker mode fallback requires a numeric mode: $mode_value" >&2
            return 1
        }
        docker_host_path_operation chmod "$mode_value" "$@" && return 0
    fi
    echo "Failed to chmod $mode_value: $*" >&2
    return 1
}

# Change host-path ownership recursively. Rootless Podman can operate on mapped
# IDs through its user namespace. Rootful Docker can perform the same change
# through a tightly scoped bind mount without requiring host sudo access.
chown_with_engine_fallback() {
    local owner="$1"
    shift

    if chown -R "$owner" "$@" 2>/dev/null; then
        return 0
    fi
    if [ "$engine" = "podman" ] && podman unshare chown -R "$owner" "$@"; then
        return 0
    fi
    if [ "$engine" = "docker" ]; then
        [[ "$owner" =~ ^[0-9]+(:[0-9]+)?$ ]] || {
            echo "Docker ownership fallback requires a numeric UID or UID:GID: $owner" >&2
            return 1
        }
        docker_host_path_operation chown "$owner" "$@" && return 0
    fi
    echo "Failed to chown $owner: $*" >&2
    return 1
}

# Apply an application-owned host permission specification consistently.
# Each argument is one `path|owner|mode` entry. Use `-` when ownership or mode
# must remain unchanged. Missing paths are skipped so optional bind sources can
# be declared alongside required ones without wrapper-owned condition loops.
apply_host_permission_spec() {
    local entry path owner mode extra

    for entry in "$@"; do
        IFS='|' read -r path owner mode extra <<< "$entry"
        if [ -z "$path" ] || [ -z "$owner" ] || [ -z "$mode" ] || [ -n "$extra" ]; then
            echo "Invalid host permission entry '$entry'; expected path|owner|mode." >&2
            return 1
        fi
        [ -e "$path" ] || continue
        if [ "$owner" != "-" ]; then
            chown_with_engine_fallback "$owner" "$path" || return 1
        fi
        if [ "$mode" != "-" ]; then
            chmod_with_engine_fallback "$mode" "$path" || return 1
        fi
    done
}

# Create and repair application-writable directories as seen inside a running
# container. Each argument is one `path|owner|mode` entry. Ownership and mode
# are applied recursively because mounted directory contents may pre-exist.
# Arguments: container ID/name, followed by directory specifications.
apply_container_directory_permission_spec() {
    local container_id="$1"
    shift
    local entry path owner mode extra

    for entry in "$@"; do
        IFS='|' read -r path owner mode extra <<< "$entry"
        if [ -z "$path" ] || [ -z "$owner" ] || [ -z "$mode" ] || [ -n "$extra" ]; then
            echo "Invalid container directory permission entry '$entry'; expected path|owner|mode." >&2
            return 1
        fi
        engine_exec exec --user 0 "$container_id" sh -c '
            path="$1"
            owner="$2"
            mode="$3"
            mkdir -p "$path"
            chown -R "$owner" "$path"
            chmod -R "$mode" "$path"
        ' _ "$path" "$owner" "$mode" || return 1
    done
}

# Copy an installation configuration into a running container and immediately
# restrict it to the application identity. The destination parent must already
# exist as part of the staged application source.
# Arguments: container ID/name, host file, container path, UID, GID.
stage_container_runtime_config() {
    local container_id="$1"
    local host_path="$2"
    local container_path="$3"
    local application_uid="$4"
    local application_gid="$5"

    if [ ! -f "$host_path" ]; then
        echo "Runtime installation configuration not found: $host_path" >&2
        return 1
    fi
    if [ -z "$container_path" ] || [ "$container_path" = "/" ]; then
        echo "Invalid runtime installation configuration destination: '$container_path'" >&2
        return 1
    fi
    engine_exec cp "$host_path" "${container_id}:$container_path" || return 1
    engine_exec exec --user 0 "$container_id" \
        chown "$application_uid:$application_gid" "$container_path" || return 1
    engine_exec exec --user 0 "$container_id" chmod 0600 "$container_path"
}

# Remove the exact temporary runtime configuration staged for bootstrap.
# Arguments: container ID/name, container path.
remove_container_runtime_config() {
    local container_id="$1"
    local container_path="$2"

    if [ -z "$container_path" ] || [ "$container_path" = "/" ]; then
        echo "Refusing to remove invalid runtime configuration path: '$container_path'" >&2
        return 1
    fi
    engine_exec exec --user 0 "$container_id" rm -f -- "$container_path"
}

# Run an application's smoke-test script with the standard deployment context.
# Arguments: script, mode, engine, Compose file, optional environment file,
# followed by application-specific smoke-test arguments.
run_standard_smoke_test() {
    local smoke_script="$1"
    local deployment_mode="$2"
    local container_engine="$3"
    local compose_path="$4"
    local environment_path="$5"
    shift 5
    local -a smoke_args=(--engine "$container_engine" --compose_file "$compose_path")

    if [ ! -f "$smoke_script" ]; then
        echo "Smoke-test script not found: $smoke_script" >&2
        return 1
    fi
    if [ "$deployment_mode" = "test" ]; then
        smoke_args+=(--test)
    elif [ -n "$environment_path" ]; then
        smoke_args+=(--env_file "$environment_path")
    fi
    echo "Running deployment smoke test: $smoke_script"
    bash "$smoke_script" "${smoke_args[@]}" "$@"
}

# Convert a URL or host:port value into a valid Apache ServerName host.
normalize_apache_server_name() {
    local value="$1"

    value="${value#http://}"
    value="${value#https://}"
    value="${value%%/*}"
    value="${value%%:*}"
    if [ -z "$value" ] || [ "$value" = "*" ]; then
        value="localhost"
    fi
    echo "$value"
}

# Escape a value before placing it in the replacement side of a sed command.
sed_replacement_escape() {
    printf '%s' "$1" | sed -e 's/[\\&|]/\\&/g'
}

# Escape a literal value before placing it in the search side of a sed command.
sed_search_escape() {
    printf '%s' "$1" | sed -e 's/[][\\.^$*+?{}()|]/\\&/g'
}

# Read one shell-style installation setting in an isolated Bash process.
# Isolation prevents configuration assignments from changing installer state.
read_install_conf_value() {
    local key="$1"
    local file="$2"
    [ -f "$file" ] || return 0
    bash -c '
        set -a
        . "$1"
        key="$2"
        printf "%s" "${!key-}"
    ' _ "$file" "$key"
}

# Return the first non-empty value among compatible setting names. This allows
# the shared renderer to bridge established RELECOV names and the standard
# scaffold names without duplicating application configuration files.
read_install_conf_first() {
    local file="$1"
    shift
    local key value
    for key in "$@"; do
        value="$(read_install_conf_value "$key" "$file")"
        if [ -n "$value" ]; then
            echo "$value"
            return 0
        fi
    done
}

# Return an environment override, then a configuration value, then a default.
# Arguments: setting name, configuration file, default value.
config_value_or_default() {
    local key="$1"
    local file="$2"
    local default_value="$3"
    local env_value="${!key:-}"
    local config_value=""

    if [ -n "$env_value" ]; then
        echo "$env_value"
        return 0
    fi
    config_value="$(read_install_conf_value "$key" "$file")"
    if [ -n "$config_value" ]; then
        echo "$config_value"
    else
        echo "$default_value"
    fi
}

# Return a required setting using the same environment-over-file precedence.
# This is used while building the generated Compose environment so a missing
# secret or host path fails before Compose changes deployment state.
config_value() {
    local key="$1"
    local file="$2"
    local value
    value="$(config_value_or_default "$key" "$file" "")"
    if [ -z "$value" ]; then
        echo "Required setting $key is missing from environment and $file" >&2
        return 1
    fi
    printf '%s\n' "$value"
}

# Normalize a bind-mount setting to a file path. Operators may configure either
# the final file or its parent directory; the supplied default is app-owned.
normalize_bind_file_path() {
    local value="$1"
    local default_path="$2"
    local filename="$3"
    local expected_suffix=".${filename##*.}"

    if [ -z "$value" ]; then
        echo "$default_path"
    elif [ -d "$value" ] || [[ "$value" = */ ]] || [[ "$value" != *"$expected_suffix" ]]; then
        echo "${value%/}/$filename"
    else
        echo "$value"
    fi
}

# Create a configuration file from a plain-text template. Arguments are:
# source template, destination file, destination mode, then repeating literal
# PLACEHOLDER VALUE pairs. Every placeholder occurrence is replaced by its
# value; values are not executed or sourced as shell code. The completed file
# is installed through the Docker/Podman-compatible host file helpers.
render_config_template() {
    local src="$1"
    local dst="$2"
    local file_mode="$3"
    shift 3
    local tmp_file token value escaped_token escaped_value

    if [ $(( $# % 2 )) -ne 0 ]; then
        echo "render_config_template requires PLACEHOLDER VALUE pairs" >&2
        return 2
    fi

    tmp_file="$(mktemp)"
    cp "$src" "$tmp_file"
    while [ "$#" -gt 0 ]; do
        token="$1"
        value="$2"
        shift 2
        escaped_token="$(sed_search_escape "$token")"
        escaped_value="$(sed_replacement_escape "$value")"
        sed -E -i "s|$escaped_token|$escaped_value|g" "$tmp_file"
    done
    if copy_with_podman_fallback "$tmp_file" "$dst" \
        && chmod_with_engine_fallback "$file_mode" "$dst"; then
        rm -f "$tmp_file"
        return 0
    fi
    rm -f "$tmp_file"
    return 1
}

# Render every ${UPPER_CASE_VARIABLE} found in a repository-owned configuration
# source from the already loaded deployment environment. This keeps application
# topology editable in normal configuration files while ensuring Compose mounts
# completed files with no unresolved deployment placeholders.
# Arguments: source configuration, destination file, destination mode.
render_environment_config_template() {
    local src="$1"
    local dst="$2"
    local file_mode="$3"
    local token variable
    local -a replacements=()

    [ -f "$src" ] || {
        echo "Configuration source not found: $src" >&2
        return 1
    }
    while IFS= read -r token; do
        [ -n "$token" ] || continue
        variable="${token#\$\{}"
        variable="${variable%\}}"
        if ! [[ -v "$variable" ]]; then
            echo "Required template variable $variable is not set for $src" >&2
            return 1
        fi
        replacements+=("$token" "${!variable}")
    done < <(grep -oE '\$\{[A-Z][A-Z0-9_]*\}' "$src" | sort -u || true)

    render_config_template "$src" "$dst" "$file_mode" "${replacements[@]}"
}

# Quote one value for literal use in a Compose environment file. Single-quoted
# dotenv values are not interpolated; embedded apostrophes are escaped.
compose_environment_quote() {
    local value="$1"
    if [[ "$value" == *$'\n'* || "$value" == *$'\r'* ]]; then
        echo "Compose environment values must not contain newlines." >&2
        return 1
    fi
    value="${value//\'/\\\'}"
    printf "'%s'" "$value"
}

# Generate one protected Compose environment file from application settings.
# Arguments: output path, settings-source array name, explicit-value array name.
#
# Settings sources contain `PREFIX|path` entries. Every uppercase assignment in
# the file is emitted as PREFIX_KEY, or as KEY when PREFIX is empty. Explicit
# values contain `KEY|value` entries for derived values such as image names and
# Git revisions. Duplicate output keys, invalid names, and multiline values are
# rejected. The result is installed atomically with mode 0600.
write_compose_environment_file() {
    local output_path="$1"
    local sources_name="$2"
    local values_name="$3"
    local -n settings_sources_ref="$sources_name"
    local -n explicit_values_ref="$values_name"
    local output_dir temporary_file entry prefix settings_path key output_key value quoted_value
    local -A emitted_keys=()

    output_dir="$(dirname "$output_path")"
    [ -d "$output_dir" ] || {
        echo "Compose environment output directory not found: $output_dir" >&2
        return 1
    }
    temporary_file="$(mktemp "$output_dir/.compose-env.XXXXXX")" || return 1
    chmod 0600 "$temporary_file" || { rm -f "$temporary_file"; return 1; }

    for entry in "${settings_sources_ref[@]}"; do
        if [[ "$entry" != *"|"* ]]; then
            echo "Invalid Compose settings source '$entry'; expected PREFIX|path." >&2
            rm -f "$temporary_file"
            return 1
        fi
        prefix="${entry%%|*}"
        settings_path="${entry#*|}"
        if [ -n "$prefix" ] && [[ ! "$prefix" =~ ^[A-Z_][A-Z0-9_]*$ ]]; then
            echo "Invalid Compose environment prefix: $prefix" >&2
            rm -f "$temporary_file"
            return 1
        fi
        if [ ! -f "$settings_path" ]; then
            echo "Compose settings source not found: $settings_path" >&2
            rm -f "$temporary_file"
            return 1
        fi

        while IFS= read -r key; do
            [ -n "$key" ] || continue
            if [ -n "$prefix" ]; then output_key="${prefix}_${key}"; else output_key="$key"; fi
            if [ -n "${emitted_keys[$output_key]:-}" ]; then
                echo "Duplicate Compose environment variable: $output_key" >&2
                rm -f "$temporary_file"
                return 1
            fi
            value="$(read_install_conf_value "$key" "$settings_path")"
            quoted_value="$(compose_environment_quote "$value")" \
                || { rm -f "$temporary_file"; return 1; }
            printf '%s=%s\n' "$output_key" "$quoted_value" >> "$temporary_file" \
                || { rm -f "$temporary_file"; return 1; }
            emitted_keys["$output_key"]=1
        done < <(sed -nE \
            's/^[[:space:]]*(export[[:space:]]+)?([A-Z_][A-Z0-9_]*)[[:space:]]*=.*/\2/p' \
            "$settings_path" | sort -u)
    done

    for entry in "${explicit_values_ref[@]}"; do
        if [[ "$entry" != *"|"* ]]; then
            echo "Invalid explicit Compose value '$entry'; expected KEY|value." >&2
            rm -f "$temporary_file"
            return 1
        fi
        output_key="${entry%%|*}"
        value="${entry#*|}"
        if [[ ! "$output_key" =~ ^[A-Z_][A-Z0-9_]*$ ]]; then
            echo "Invalid Compose environment variable name: $output_key" >&2
            rm -f "$temporary_file"
            return 1
        fi
        if [ -n "${emitted_keys[$output_key]:-}" ]; then
            echo "Duplicate Compose environment variable: $output_key" >&2
            rm -f "$temporary_file"
            return 1
        fi
        quoted_value="$(compose_environment_quote "$value")" \
            || { rm -f "$temporary_file"; return 1; }
        printf '%s=%s\n' "$output_key" "$quoted_value" >> "$temporary_file" \
            || { rm -f "$temporary_file"; return 1; }
        emitted_keys["$output_key"]=1
    done

    mv -f "$temporary_file" "$output_path" || { rm -f "$temporary_file"; return 1; }
    chmod 0600 "$output_path"
}

# Load the protected dotenv file generated by write_compose_environment_file.
# Direct image builds and host preparation then consume exactly the same
# prefixed deployment values that Compose interpolates later.
load_compose_environment_file() {
    local file="$1"
    local disable_allexport_after="true"
    [ -f "$file" ] || {
        echo "Compose environment file not found: $file" >&2
        return 1
    }
    [[ $- == *a* ]] && disable_allexport_after="false"
    set -a
    # shellcheck disable=SC1090
    source "$file"
    [ "$disable_allexport_after" = false ] || set +a
}

# Validate an installation configuration and, when it is outside the selected
# build context, copy it to a temporary context file. All results are explicit
# named output variables; applications still choose service contexts/defaults.
# Arguments: input path, base directory, build context, service, mode,
#            host-path output name, context-relative output name,
#            temporary-path output name.
prepare_install_configuration() {
    local input_path="$1"
    local base_directory="$2"
    local build_context="$3"
    local service_name="$4"
    local deployment_mode="$5"
    local -n host_path_result="$6"
    local -n context_path_result="$7"
    local -n temporary_path_result="$8"
    local resolved_host_path=""
    local resolved_context=""
    local temporary_name=""

    host_path_result=""
    context_path_result=""
    temporary_path_result=""

    if [ ! -d "$build_context" ]; then
        echo "Build context directory '$build_context' for service '$service_name' not found" >&2
        return 1
    fi
    resolved_context="$(cd "$build_context" && pwd -P)"

    if [[ "$input_path" = /* ]]; then
        resolved_host_path="$input_path"
    else
        resolved_host_path="${base_directory%/}/$input_path"
    fi
    if [ ! -f "$resolved_host_path" ]; then
        echo "Install configuration '$input_path' for service '$service_name' not found" >&2
        return 1
    fi
    resolved_host_path="$(cd "$(dirname "$resolved_host_path")" && pwd -P)/$(basename "$resolved_host_path")"

    if [ "$deployment_mode" = "production" ] \
        && [[ "$(basename "$resolved_host_path")" != *settings*.txt ]]; then
        echo "Production configuration filenames must match *settings*.txt so .dockerignore excludes them from COPY: $resolved_host_path" >&2
        return 1
    fi

    if [[ "$resolved_host_path" != "$resolved_context/"* ]]; then
        if [ "$deployment_mode" = "production" ]; then
            temporary_name=".tmp_docker_install_conf_${service_name}_$$.txt"
        else
            # Test configuration is explicitly non-sensitive and must remain in
            # the context because test builds do not receive a build secret.
            temporary_name=".tmp_docker_test_install_conf_${service_name}_$$.txt"
        fi
        temporary_path_result="$resolved_context/$temporary_name"
        echo "Copying $resolved_host_path into temporary file $temporary_path_result for service '$service_name'." >&2
        cp "$resolved_host_path" "$temporary_path_result" || return 1
        resolved_host_path="$temporary_path_result"
    fi

    host_path_result="$resolved_host_path"
    context_path_result="${resolved_host_path#$resolved_context/}"
}

# Remove explicitly supplied temporary files. It is intended for EXIT traps;
# callers retain ownership of the list and no directory is ever removed.
cleanup_files() {
    local file
    for file in "$@"; do
        if [ -n "$file" ] && [ -f "$file" ]; then
            rm -f "$file"
        fi
    done
}

# Return success when the first argument exactly matches one of the remaining
# arguments. An empty value list is valid and returns failure.
# Arguments: value to find, zero or more candidate values.
array_contains() {
    local wanted="$1"
    shift
    local candidate
    for candidate in "$@"; do
        if [ "$candidate" = "$wanted" ]; then
            return 0
        fi
    done
    return 1
}

# Test whether a service is declared in the selected Compose model. A stopped
# service still exists, so this deliberately uses config rather than ps.
service_exists() {
    compose_with_env_exec -f "$compose_file" config --services 2>/dev/null \
        | grep -Fxq "$1"
}

# Resolve a service to a container name/ID, supporting application legacy names
# first. Otherwise inspect every container returned by the selected Compose
# project and match its service label. Listing the project first avoids both a
# cross-project "app" collision and the unsupported `ps -q SERVICE` syntax in
# podman-compose 1.0.x.
resolve_service_container() {
    local service_name="$1"
    local service_container=""
    local container_name=""

    # Applications may define service_container_name to support explicit legacy
    # container_name values. Compose labels are the generic fallback.
    if declare -F service_container_name >/dev/null 2>&1; then
        container_name="$(service_container_name "$service_name")"
    fi
    if [ -n "$container_name" ] \
        && engine_exec inspect -f '{{.Id}}' "$container_name" >/dev/null 2>&1; then
        service_container="$container_name"
    elif [ -n "${compose_file:-}" ]; then
        local candidate candidate_service
        while IFS= read -r candidate; do
            [ -n "$candidate" ] || continue
            candidate_service="$(engine_exec inspect -f \
                '{{ index .Config.Labels "com.docker.compose.service" }}' \
                "$candidate" 2>/dev/null || true)"
            if [ "$candidate_service" = "$service_name" ]; then
                service_container="$candidate"
                break
            fi
        done < <(compose_with_env_exec -f "$compose_file" ps -q 2>/dev/null || true)
    fi
    # Retain the label query as a compatibility fallback for callers that do
    # not have a Compose file in scope (including older application wrappers).
    if [ -z "$service_container" ]; then
        service_container="$(engine_exec ps -a \
            --filter "label=com.docker.compose.service=${service_name}" \
            --format '{{.ID}}' | head -n 1)"
    fi
    if [ -z "$service_container" ]; then
        echo "Error: unable to resolve container ID for service '$service_name'." >&2
        return 1
    fi
    echo "$service_container"
}

# Fail with useful logs unless a resolved service container is running.
ensure_service_running() {
    local service_name="$1"
    local service_container="$2"

    if ! engine_exec inspect -f '{{.State.Running}}' "$service_container" >/dev/null 2>&1; then
        echo "Error: service '$service_name' container does not exist." >&2
        exit 1
    fi
    if [ "$(engine_exec inspect -f '{{.State.Running}}' "$service_container")" != "true" ]; then
        echo "Error: service '$service_name' container is not running. Showing logs:" >&2
        engine_exec logs --tail 200 "$service_container" >&2 || true
        exit 1
    fi
    echo "$service_container"
}
