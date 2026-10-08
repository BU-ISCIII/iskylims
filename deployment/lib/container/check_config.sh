#!/usr/bin/env bash
# Cross-check generated deployment wiring before images are built.
#
# This file is centrally managed. Applications must not edit vendored copies.
#
# `container_install.sh` calls this checker after Compose validation. It reads
# only generated artifacts: the protected Compose environment file, the Compose
# file and, when the Apache add-on is selected, the rendered Apache
# configuration. Each finding names a rule and the setting to change. Values
# are printed only for URL and host settings, never for other settings.
#
# Usage: check_config.sh --mode {production,test} --env-file FILE
#            --compose-file FILE [--apache-config-dir DIR]
#            [--service NAME=PROFILE]...
#
# Rules:
#   unresolved-placeholder  a CHANGE_ME marker remains in a setting
#   unknown-host            an internal URL or DB_HOST names no Compose service
#   upstream-port           an internal URL port differs from the target's port
#   loopback-target         a server-side proxy target points at the container itself
#   database-host           DB_HOST ignores the Compose-managed database service
#   django-hostname         a host name with underscores reaches Django
#   allowed-hosts           a host reaching Django is missing from ALLOWED_HOSTS
#   keycloak-url            an OIDC issuer or JWKS URL is malformed
#   keycloak-origin         Keycloak public URLs disagree across settings
#   keycloak-realm          Keycloak realms disagree across settings
#   keycloak-vhost          Keycloak public URL is not routed to Keycloak by Apache
#   canonical-url           AUTH_URL and NEXTAUTH_URL disagree
#
# Production findings are errors and fail the installation. Test findings are
# reported as warnings because test settings intentionally use shortcuts.
# Findings the checker cannot confirm (for example a single-label host that may
# resolve through institutional DNS) are always warnings.
#
# The checker needs only Bash 4.4 or newer and `sort`, so target hosts do not
# need any other interpreter.

set -euo pipefail
export LC_ALL=C
shopt -s nullglob

if ((BASH_VERSINFO[0] < 4 || (BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 4))); then
    echo "Error: check_config.sh requires Bash 4.4 or newer" >&2
    exit 2
fi

readonly PLACEHOLDER="CHANGE_ME"
readonly URL_PATTERN="https?://[^[:space:],;'\"<>]+"
readonly ISSUER_PATH='^(/auth)?/realms/([^/]+)/?$'
readonly JWKS_PATH='^(/auth)?/realms/([^/]+)/protocol/openid-connect/certs/?$'
readonly NL=$'\n' TAB=$'\t'
KEYCLOAK_ORIGIN_SUFFIXES=(KEYCLOAK_URL KEYCLOAK_PUBLIC_URL KEYCLOAK_ADMIN_API_BASE_URL)
KEYCLOAK_REALM_SUFFIXES=(KEYCLOAK_REALM KEYCLOAK_ADMIN_API_REALM)
# Server-side Keycloak calls may use an internal route instead of the public URL.
KEYCLOAK_SERVER_SIDE_SUFFIXES=(OIDC_JWKS_URL KEYCLOAK_ADMIN_API_BASE_URL)

# Settings from the Compose environment file; ENV_KEYS holds them sorted.
declare -A ENV=()
ENV_KEYS=()
# Compose services in file order. Names and ports are newline-delimited sets.
SERVICES=()
declare -A SERVICE_NAMES=() SERVICE_PORTS=()
# Application services and their profiles in command-line order.
declare -A PROFILES=()
PROFILE_ORDER=()
# Compose file currently being validated.
COMPOSE_FILE=""
# Apache ProxyPass targets with their virtual host names (newline-delimited).
ROUTE_SOURCES=() ROUTE_NAMES=() ROUTE_PRESERVE=() ROUTE_TARGETS=()
# Host headers each Django service receives, as "host<TAB>source" lines.
declare -A ARRIVALS=()
FINDING_RULES=() FINDING_MESSAGES=() FINDING_CONFIRMED=()

# Helpers return their result in REPLY to avoid subshells.

print_help() {
    local line
    while IFS= read -r line; do
        case "$line" in
            '#!'*) continue ;;
            '#'*) line="${line#\#}"; printf '%s\n' "${line# }" ;;
            *) break ;;
        esac
    done < "${BASH_SOURCE[0]}"
}

usage_error() {
    echo "usage: check_config.sh --mode {production,test} --env-file FILE" \
        "--compose-file FILE [--apache-config-dir DIR] [--service NAME=PROFILE]..." >&2
    echo "check_config.sh: error: $1" >&2
    exit 2
}

add_finding() {
    FINDING_RULES+=("$1")
    FINDING_MESSAGES+=("$2")
    FINDING_CONFIRMED+=("${3:-1}")
}

trim() {
    local value="$1"
    value="${value#"${value%%[![:space:]]*}"}"
    REPLY="${value%"${value##*[![:space:]]}"}"
}

env_prefix() {
    REPLY="${1^^}"
    REPLY="${REPLY//[^A-Z0-9]/_}"
}

normalized_name() {
    REPLY="${1,,}"
    REPLY="${REPLY//[_.]/-}"
}

is_loopback() {
    case "$1" in localhost|127.0.0.1|::1) return 0 ;; esac
    return 1
}

ends_with_any() {
    local value="$1" suffix
    shift
    for suffix in "$@"; do
        if [[ "$value" == *"$suffix" ]]; then return 0; fi
    done
    return 1
}

# Quote a value the way Python's repr() does, as earlier checker output did.
py_repr() {
    local value="$1" backslash='\' quote="'"
    value="${value//"$backslash"/"$backslash$backslash"}"
    value="${value//"$TAB"/"\\t"}"
    if [[ "$value" == *"$quote"* && "$value" != *'"'* ]]; then
        REPLY="\"$value\""
    else
        REPLY="'${value//"$quote"/"$backslash$quote"}'"
    fi
}

# Recognize an IPv6 address. Callers already treat any host containing a dot
# (IPv4 addresses included) as external, so only dotless IPv6 forms matter.
is_ipv6() {
    local host="${1%%\%*}" group
    local -a groups
    [[ "$host" == *:* && "$host" =~ ^[0-9A-Fa-f:]+$ ]] || return 1
    [[ "$host" != *::*::* ]] || return 1
    [[ "$host" != :* || "$host" == ::* ]] || return 1
    [[ "$host" != *: || "$host" == *:: ]] || return 1
    IFS=: read -ra groups <<< "${host/::/:}"
    for group in "${groups[@]}"; do
        [ "${#group}" -le 4 ] || return 1
    done
    if [[ "$host" == *::* ]]; then
        [ "${#groups[@]}" -le 8 ]
    else
        [ "${#groups[@]}" -eq 8 ]
    fi
}

default_port() {
    case "$1" in
        http) REPLY=80 ;;
        https) REPLY=443 ;;
        *) REPLY="" ;;
    esac
}

# Split a URL like urllib.parse.urlsplit. Sets U_SCHEME, U_NETLOC, U_HOST
# (lowercase), U_PORT (empty when absent), U_PORT_VALID (0 for a malformed or
# out-of-range port) and U_PATH. Callers that must keep their own split
# declare these names local before calling.
url_split() {
    local rest="$1" hostinfo port="" open="[" close="]"
    U_SCHEME="" U_NETLOC="" U_HOST="" U_PORT="" U_PORT_VALID=1 U_PATH=""
    if [[ "$rest" =~ ^([A-Za-z][A-Za-z0-9+.-]*):(.*)$ ]]; then
        U_SCHEME="${BASH_REMATCH[1],,}"
        rest="${BASH_REMATCH[2]}"
    fi
    if [[ "$rest" == //* ]]; then
        rest="${rest#//}"
        U_NETLOC="${rest%%[/?#]*}"
        rest="${rest:${#U_NETLOC}}"
    fi
    U_PATH="${rest%%[?#]*}"

    hostinfo="${U_NETLOC##*@}"
    if [[ "$hostinfo" == *"$open"* ]]; then
        hostinfo="${hostinfo#*"$open"}"
        U_HOST="${hostinfo%%"$close"*}"
        if [[ "$hostinfo" == *"$close"*:* ]]; then port="${hostinfo#*"$close"}"; port="${port#*:}"; fi
    else
        U_HOST="${hostinfo%%:*}"
        if [[ "$hostinfo" == *:* ]]; then port="${hostinfo#*:}"; fi
    fi
    U_HOST="${U_HOST,,}"

    if [ -n "$port" ]; then
        if [[ "$port" =~ ^[0-9]+$ ]]; then
            while [[ "$port" == 0?* ]]; do port="${port#0}"; done
            if [ "${#port}" -le 5 ] && [ "$port" -le 65535 ]; then
                U_PORT="$port"
            else
                U_PORT_VALID=0
            fi
        else
            U_PORT_VALID=0
        fi
    fi
}

# Remove credentials before a URL is printed.
display_url() {
    local U_SCHEME U_NETLOC U_HOST U_PORT U_PORT_VALID U_PATH
    url_split "$1"
    REPLY="$1"
    if [[ "$U_NETLOC" == *@* ]]; then
        REPLY="${1/"$U_NETLOC"/"${U_NETLOC##*@}"}"
    fi
}

# Return the bare host of an Apache ServerName or Host-like value.
host_name() {
    local value colons
    trim "$1"
    value="${REPLY,,}"
    if [[ "$value" == *://* ]]; then
        local U_SCHEME U_NETLOC U_HOST U_PORT U_PORT_VALID U_PATH
        url_split "$value"
        value="$U_HOST"
    fi
    colons="${value//[^:]/}"
    if [[ "$value" != "["* && "${#colons}" -eq 1 ]]; then
        value="${value%:*}"
    fi
    REPLY="$value"
}

origin() {
    local U_SCHEME U_NETLOC U_HOST U_PORT U_PORT_VALID U_PATH
    url_split "$1"
    [ "$U_PORT_VALID" = 1 ] || U_PORT=""
    default_port "$U_SCHEME"
    if [ -z "$U_PORT" ] || [ "$U_PORT" = "$REPLY" ]; then
        REPLY="$U_SCHEME://$U_HOST"
    else
        REPLY="$U_SCHEME://$U_HOST:$U_PORT"
    fi
}

# Mirror django.http.request.validate_host for ALLOWED_HOSTS matching.
# Arguments: host, then the ALLOWED_HOSTS patterns.
django_host_allowed() {
    local host="$1" pattern
    shift
    while [[ "$host" == *. ]]; do host="${host%.}"; done
    host="${host,,}"
    for pattern in "$@"; do
        pattern="${pattern,,}"
        if [ "$pattern" = "*" ] || [ "$pattern" = "$host" ]; then return 0; fi
        if [[ "$pattern" == .* && ( "$host" == *"$pattern" || "$host" == "${pattern#.}" ) ]]; then
            return 0
        fi
    done
    return 1
}

# Read the dotenv file written by write_compose_environment_file.
read_environment() {
    local line key raw
    while IFS= read -r line || [ -n "$line" ]; do
        line="${line%$'\r'}"
        trim "$line"
        if [ -z "$REPLY" ] || [[ "$REPLY" == "#"* ]] || [[ "$line" != *=* ]]; then
            continue
        fi
        trim "${line%%=*}"
        key="$REPLY"
        raw="${line#*=}"
        if [ "${#raw}" -ge 2 ] && [[ "$raw" == \'*\' ]]; then
            raw="${raw:1:${#raw}-2}"
            raw="${raw//"\\'"/"'"}"
        fi
        [ -z "$key" ] || ENV["$key"]="$raw"
    done < "$1"
    if [ "${#ENV[@]}" -gt 0 ]; then
        mapfile -t ENV_KEYS < <(printf '%s\n' "${!ENV[@]}" | sort)
    fi
}

yaml_scalar() {
    local value="${1%% #*}"
    trim "$value"
    value="$REPLY"
    while [[ "$value" == [\"\']* ]]; do value="${value:1}"; done
    while [[ "$value" == *[\"\'] ]]; do value="${value%?}"; done
    REPLY="$value"
}

# Arguments: service, list key (`expose`, `aliases` or `container_name`), value.
add_compose_value() {
    local service="$1" key="$2" value="$3"
    if [ -z "$value" ] || [[ "$value" == *'$'* ]]; then return 0; fi
    if [ "$key" = "expose" ]; then
        if [[ "$value" =~ ^([0-9]+) ]]; then
            value="$((10#${BASH_REMATCH[1]}))"
            if [[ "${SERVICE_PORTS[$service]}" != *"$NL$value$NL"* ]]; then
                SERVICE_PORTS["$service"]+="$value$NL"
            fi
        fi
    elif [[ "${SERVICE_NAMES[$service]}" != *"$NL${value,,}$NL"* ]]; then
        SERVICE_NAMES["$service"]+="${value,,}$NL"
    fi
}

# Read service names, network aliases, container names and exposed ports.
# Generated Compose files use conventional block indentation, so a small line
# scanner is sufficient and keeps the checker free of YAML tools.
read_compose_services() {
    local raw text indent in_services=0 current="" list_key="" key rest inner item
    local quote="[\"']"
    local service_pattern="^${quote}?([A-Za-z0-9._-]+)${quote}?:\$"
    local key_pattern='^([A-Za-z0-9_]+):[[:space:]]*(.*)$'
    local -a items
    while IFS= read -r raw || [ -n "$raw" ]; do
        trim "$raw"
        text="$REPLY"
        if [ -z "$text" ] || [[ "$text" == "#"* ]]; then continue; fi
        inner="${raw%%[! ]*}"
        indent="${#inner}"
        if [ "$indent" -eq 0 ]; then
            if [ "$text" = "services:" ]; then in_services=1; else in_services=0; fi
            current=""
            continue
        fi
        [ "$in_services" -eq 1 ] || continue
        if [ "$indent" -eq 2 ] && [[ "$text" =~ $service_pattern ]]; then
            current="${BASH_REMATCH[1]}"
            if [ -z "${SERVICE_NAMES[$current]+set}" ]; then
                SERVICES+=("$current")
                SERVICE_NAMES["$current"]="$NL${current,,}$NL"
                SERVICE_PORTS["$current"]="$NL"
            fi
            list_key=""
            continue
        fi
        [ -n "$current" ] || continue
        if [[ "$text" == -* ]]; then
            if [ -n "$list_key" ]; then
                yaml_scalar "${text:1}"
                add_compose_value "$current" "$list_key" "$REPLY"
            fi
            continue
        fi
        [[ "$text" =~ $key_pattern ]] || continue
        key="${BASH_REMATCH[1]}"
        rest="${BASH_REMATCH[2]}"
        list_key=""
        if [ "$key" = "container_name" ]; then
            yaml_scalar "$rest"
            add_compose_value "$current" "$key" "$REPLY"
        elif [ "$key" = "aliases" ] || [ "$key" = "expose" ]; then
            if [[ "$rest" == "["* ]]; then
                inner="$rest"
                while [[ "$inner" == [][]* ]]; do inner="${inner:1}"; done
                while [[ "$inner" == *[][] ]]; do inner="${inner%?}"; done
                IFS=, read -ra items <<< "$inner"
                for item in "${items[@]}"; do
                    yaml_scalar "$item"
                    add_compose_value "$current" "$key" "$REPLY"
                done
            elif [ -z "$rest" ]; then
                list_key="$key"
            fi
        fi
    done < "$1"
}

# Read ProxyPass targets and the host names forwarded to them.
read_apache_routes() {
    local conf raw text number directive server_names preserve index word
    local -a words pending_lines pending_targets
    for conf in "$1"/*.conf; do
        [ -f "$conf" ] || continue
        server_names="" preserve=0 number=0
        pending_lines=() pending_targets=()
        while IFS= read -r raw || [ -n "$raw" ]; do
            number=$((number + 1))
            trim "$raw"
            text="$REPLY"
            if [ -z "$text" ] || [[ "$text" == "#"* ]]; then continue; fi
            read -ra words <<< "$text"
            directive="${words[0],,}"
            if [[ "$directive" == "<virtualhost"* ]]; then
                server_names="" preserve=0
                pending_lines=() pending_targets=()
            elif [ "$directive" = "servername" ] || [ "$directive" = "serveralias" ]; then
                for word in "${words[@]:1}"; do server_names+="$word$NL"; done
            elif [ "$directive" = "proxypreservehost" ] && [ "${#words[@]}" -gt 1 ]; then
                if [ "${words[1],,}" = "on" ]; then preserve=1; else preserve=0; fi
            elif [ "$directive" = "proxypass" ] && [ "${#words[@]}" -gt 2 ] \
                    && [ "${words[2]}" != "!" ]; then
                pending_lines+=("$number")
                pending_targets+=("${words[2]}")
            elif [ "$directive" = "</virtualhost>" ]; then
                for index in "${!pending_targets[@]}"; do
                    ROUTE_SOURCES+=("${conf##*/}:${pending_lines[$index]}")
                    ROUTE_NAMES+=("$server_names")
                    ROUTE_PRESERVE+=("$preserve")
                    ROUTE_TARGETS+=("${pending_targets[$index]}")
                done
                pending_lines=() pending_targets=()
            fi
        done < "$conf"
    done
}

# Resolve a host to the Compose service that owns it (REPLY; empty if none).
resolve_service() {
    local host="${1,,}" service
    for service in "${SERVICES[@]}"; do
        if [[ "${SERVICE_NAMES[$service]}" == *"$NL$host$NL"* ]]; then
            REPLY="$service"
            return 0
        fi
    done
    REPLY=""
    return 1
}

# Name the Compose service whose name differs from host only in `_`/`.`/`-`.
suggest_service() {
    local wanted service name
    local -a names
    normalized_name "$1"
    wanted="$REPLY"
    for service in "${SERVICES[@]}"; do
        mapfile -t names <<< "${SERVICE_NAMES[$service]}"
        for name in "${names[@]}"; do
            [ -n "$name" ] || continue
            normalized_name "$name"
            if [ "$REPLY" = "$wanted" ]; then
                REPLY="$service"
                return 0
            fi
        done
    done
    REPLY=""
}

is_internal() {
    local host="${1,,}"
    is_loopback "$host" || [ "$host" = "host.docker.internal" ] || resolve_service "$host"
}

is_django() {
    [ "${PROFILES[$1]:-}" = "django" ]
}

# The settings-owned APP_PORT is the service's only listener.
apply_application_ports() {
    local name port
    for name in "${PROFILE_ORDER[@]}"; do
        env_prefix "$name"
        port="${ENV[${REPLY}_APP_PORT]:-}"
        if [ -n "${SERVICE_NAMES[$name]+set}" ] && [[ "$port" =~ ^[0-9]+$ ]]; then
            SERVICE_PORTS["$name"]="$NL$((10#$port))$NL"
        fi
    done
}

# Resolve an internal URL host and validate its port. The caller has already
# split the URL with url_split. Sets FOUND_SERVICE to the resolved service.
# Arguments: label, URL, lowercase host.
check_internal_host() {
    local label="$1" url="$2" host="$3" shown host_shown suggestion port port_target expected port_value
    FOUND_SERVICE=""
    display_url "$url"
    py_repr "$REPLY"
    shown="$REPLY"
    py_repr "$host"
    host_shown="$REPLY"
    if resolve_service "$host"; then
        FOUND_SERVICE="$REPLY"
        port_target="$REPLY"
    else
        if [[ "$host" == *.* ]] || is_loopback "$host" || is_ipv6 "$host"; then
            return 0
        fi
        suggest_service "$host"
        suggestion="$REPLY"
        if [ -z "$suggestion" ]; then
            add_finding unknown-host "$label=$shown: host $host_shown is not a Compose service; confirm it resolves from inside the containers" 0
            return 0
        fi
        py_repr "$suggestion"
        add_finding unknown-host "$label=$shown: host $host_shown is not a Compose service; did you mean $REPLY?"
        # Also report a wrong port now so both fixes happen in one run.
        port_target="$suggestion"
    fi
    [ "$U_PORT_VALID" = 1 ] || return 0
    port="$U_PORT"
    if [ -z "$port" ] || [ "$port" = 0 ]; then
        default_port "$U_SCHEME"
        port="$REPLY"
    fi
    if [ "${SERVICE_PORTS[$port_target]}" != "$NL" ] \
            && [[ "${SERVICE_PORTS[$port_target]}" != *"$NL$port$NL"* ]]; then
        expected=""
        while IFS= read -r port_value; do
            [ -z "$port_value" ] || expected+="${expected:+, }$port_value"
        done < <(sort -n <<< "${SERVICE_PORTS[$port_target]}")
        add_finding upstream-port "$label=$shown: $port_target listens on $expected, not ${port:-None}"
    fi
}

check_placeholders() {
    local key
    for key in "${ENV_KEYS[@]}"; do
        if [[ "${ENV[$key]}" == *"$PLACEHOLDER"* ]]; then
            add_finding unresolved-placeholder "$key still contains $PLACEHOLDER"
        fi
    done
}

# Validate every internal URL and record the Host each Django service receives.
check_internal_urls() {
    local key value rest url host label index name
    local U_SCHEME U_NETLOC U_HOST U_PORT U_PORT_VALID U_PATH
    local -a names
    for key in "${ENV_KEYS[@]}"; do
        value="${ENV[$key]}"
        [[ "$value" != *"$PLACEHOLDER"* ]] || continue
        rest="$value"
        while [[ "$rest" =~ $URL_PATTERN ]]; do
            url="${BASH_REMATCH[0]}"
            rest="${rest#*"$url"}"
            url_split "$url"
            host="$U_HOST"
            if [ -z "$host" ] || [[ "$url" == *'$'* ]]; then continue; fi
            if [[ "$key" == *PROXY_TARGET ]] && is_loopback "$host"; then
                display_url "$url"
                py_repr "$REPLY"
                add_finding loopback-target "$key=$REPLY: inside a container $host is the calling container itself; use the Compose service name"
                continue
            fi
            check_internal_host "$key" "$url" "$host"
            if [ -n "$FOUND_SERVICE" ] && is_django "$FOUND_SERVICE"; then
                ARRIVALS["$FOUND_SERVICE"]+="$host$TAB$key$NL"
            fi
        done
    done

    for index in "${!ROUTE_TARGETS[@]}"; do
        url="${ROUTE_TARGETS[$index]}"
        url_split "$url"
        host="$U_HOST"
        if [ -z "$host" ] || [[ "$url" == *'$'* ]]; then continue; fi
        label="Apache ProxyPass ${ROUTE_SOURCES[$index]}"
        check_internal_host "$label" "$url" "$host"
        if [ -z "$FOUND_SERVICE" ] || ! is_django "$FOUND_SERVICE"; then continue; fi
        if [ "${ROUTE_PRESERVE[$index]}" = 1 ]; then
            mapfile -t names <<< "${ROUTE_NAMES[$index]}"
            for name in "${names[@]}"; do
                [ -n "$name" ] || continue
                host_name "$name"
                ARRIVALS["$FOUND_SERVICE"]+="$REPLY$TAB$label ServerName (ProxyPreserveHost On)$NL"
            done
        else
            ARRIVALS["$FOUND_SERVICE"]+="$host$TAB$label$NL"
        fi
    done
}

check_django_hosts() {
    local service_name allowed_key allowed_value has_allowed arrival host source shown item
    local -a services arrivals patterns items
    local -A seen
    [ "${#ARRIVALS[@]}" -gt 0 ] || return 0
    mapfile -t services < <(printf '%s\n' "${!ARRIVALS[@]}" | sort)
    for service_name in "${services[@]}"; do
        env_prefix "$service_name"
        allowed_key="${REPLY}_DJANGO_ALLOWED_HOSTS"
        has_allowed=0 allowed_value="" patterns=()
        if [ -n "${ENV[$allowed_key]+set}" ]; then
            has_allowed=1
            allowed_value="${ENV[$allowed_key]}"
            IFS=, read -ra items <<< "$allowed_value"
            for item in "${items[@]}"; do
                trim "$item"
                [ -z "$REPLY" ] || patterns+=("$REPLY")
            done
        fi
        seen=()
        mapfile -t arrivals <<< "${ARRIVALS[$service_name]}"
        for arrival in "${arrivals[@]}"; do
            [ -n "$arrival" ] && [ -z "${seen[$arrival]+set}" ] || continue
            seen["$arrival"]=1
            host="${arrival%%"$TAB"*}"
            source="${arrival#*"$TAB"}"
            if [ -z "$host" ] || [[ "$host" == *[*\$]* ]]; then continue; fi
            py_repr "$host"
            shown="$REPLY"
            if [[ "$host" == *_* ]]; then
                add_finding django-hostname "$source: Django rejects the Host $shown sent to $service_name because it contains an underscore; use a hyphenated service name or network alias"
            elif [ "$has_allowed" = 1 ] && [[ "$allowed_value" != *"$PLACEHOLDER"* ]] \
                    && ! django_host_allowed "$host" "${patterns[@]}"; then
                add_finding allowed-hosts "$source: $service_name receives Host $shown, which $allowed_key does not allow"
            fi
        done
    done
}


check_databases() {
    local service_name key host wanted name resolved managed first_managed host_shown suggestion
    local -a sorted_services
    mapfile -t sorted_services < <(printf '%s\n' "${SERVICES[@]}" | sort)
    for service_name in "${PROFILE_ORDER[@]}"; do
        [ "${PROFILES[$service_name]}" = "django" ] || continue
        env_prefix "$service_name"
        key="${REPLY}_DB_HOST"
        trim "${ENV[$key]:-}"
        if [ -z "$REPLY" ] || [[ "$REPLY" == *"$PLACEHOLDER"* ]]; then continue; fi
        host="${REPLY,,}"
        py_repr "$host"
        host_shown="$REPLY"
        normalized_name "$service_name"
        wanted="$REPLY-db"
        managed="$NL" first_managed=""
        for name in "${sorted_services[@]}"; do
            normalized_name "$name"
            if [ -n "$name" ] && [ "$REPLY" = "$wanted" ]; then
                managed+="$name$NL"
                first_managed="${first_managed:-$name}"
            fi
        done
        resolve_service "$host" || true
        resolved="$REPLY"
        if is_loopback "$host"; then
            add_finding database-host "$key=$host_shown: inside a container $host is the application container itself"
        elif [ -n "$first_managed" ] && [[ -z "$resolved" || "$managed" != *"$NL$resolved$NL"* ]]; then
            py_repr "$first_managed"
            add_finding database-host "$key=$host_shown: Compose manages the database service $REPLY for $service_name; use DB_HOST=$REPLY or remove that service for an external database"
        elif [ -z "$resolved" ] && [[ "$host" != *.* ]] && ! is_ipv6 "$host"; then
            suggest_service "$host"
            suggestion="$REPLY"
            if [ -n "$suggestion" ]; then
                py_repr "$suggestion"
                add_finding unknown-host "$key=$host_shown is not a Compose service; did you mean $REPLY?"
            else
                add_finding unknown-host "$key=$host_shown is not a Compose service; confirm it resolves from inside the containers" 0
            fi
        fi
    done
}

# Add a setting to the comma-separated keys recorded for a value.
# Arguments: associative array name, value, key.
record_value() {
    local -n values="$1"
    values["$2"]+="${values[$2]:+, }$3"
}

# Print "'value' in KEY, KEY; ..." for every value, sorted by value.
disagreement() {
    local -n values="$1"
    local value text=""
    while IFS= read -r value; do
        py_repr "$value"
        text+="${text:+; }$REPLY in ${values[$value]}"
    done < <(printf '%s\n' "${!values[@]}" | sort)
    REPLY="$text"
}

# Do not validate application OIDC defaults when Compose overrides them.
compose_overrides_oidc() {
    local key="$1" prefix service line in_service=0 issuer=0 jwks=0

    case "$key" in
        *_OIDC_ISSUER) prefix="${key%_OIDC_ISSUER}" ;;
        *_OIDC_JWKS_URL) prefix="${key%_OIDC_JWKS_URL}" ;;
        *) return 1 ;;
    esac
    service="${prefix,,}"
    service="${service//_/-}"

    while IFS= read -r line || [ -n "$line" ]; do
        if [[ "$line" == "  $service:" ]]; then
            in_service=1
            continue
        fi
        if [[ "$line" =~ ^\ \ [A-Za-z0-9][A-Za-z0-9_-]*: ]]; then
            in_service=0
            continue
        fi
        [ "$in_service" = 1 ] || continue
        [[ "$line" == *"OIDC_ISSUER:"*"KEYCLOAK_PUBLIC_URL"* ]] && issuer=1
        [[ "$line" == *"OIDC_JWKS_URL:"*"http://keycloak:8080/realms/"* ]] && jwks=1
    done < "$COMPOSE_FILE"

    [ "$issuer" = 1 ] && [ "$jwks" = 1 ]
}

check_keycloak() {
    local key value issuer jwks match expected server_side
    local U_SCHEME U_NETLOC U_HOST U_PORT U_PORT_VALID U_PATH
    local -A origins=() realms=()
    for key in "${ENV_KEYS[@]}"; do
        trim "${ENV[$key]}"
        value="$REPLY"
        if [ -z "$value" ] || [[ "$value" == *"$PLACEHOLDER"* ]]; then continue; fi
        issuer=0 jwks=0
        [[ "$key" != *OIDC_ISSUER ]] || issuer=1
        [[ "$key" != *OIDC_JWKS_URL ]] || jwks=1
        if { [ "$issuer" = 1 ] || [ "$jwks" = 1 ]; } \
                && compose_overrides_oidc "$key"; then
            continue
        fi
        if [ "$issuer" = 1 ] || [ "$jwks" = 1 ]; then
            url_split "$value"
            match=""
            if [ "$issuer" = 1 ] && [[ "$U_PATH" =~ $ISSUER_PATH ]]; then
                match="${BASH_REMATCH[2]}"
            elif [ "$jwks" = 1 ] && [[ "$U_PATH" =~ $JWKS_PATH ]]; then
                match="${BASH_REMATCH[2]}"
            fi
            if [ -z "$U_HOST" ] || [ -z "$match" ]; then
                if [ "$issuer" = 1 ]; then
                    expected="<keycloak-url>/realms/<realm>"
                else
                    expected="<keycloak-url>/realms/<realm>/protocol/openid-connect/certs"
                fi
                display_url "$value"
                py_repr "$REPLY"
                add_finding keycloak-url "$key=$REPLY does not match $expected"
                continue
            fi
            record_value realms "$match" "$key"
            if [ "$jwks" = 0 ] || ! is_internal "$U_HOST"; then
                origin "$value"
                record_value origins "$REPLY" "$key"
            fi
        elif ends_with_any "$key" "${KEYCLOAK_ORIGIN_SUFFIXES[@]}"; then
            url_split "$value"
            server_side=0
            ! ends_with_any "$key" "${KEYCLOAK_SERVER_SIDE_SUFFIXES[@]}" || server_side=1
            if [ -n "$U_HOST" ] && { [ "$server_side" = 0 ] || ! is_internal "$U_HOST"; }; then
                origin "$value"
                record_value origins "$REPLY" "$key"
            fi
        elif ends_with_any "$key" "${KEYCLOAK_REALM_SUFFIXES[@]}"; then
            record_value realms "$value" "$key"
        fi
    done

    if [ "${#origins[@]}" -gt 1 ]; then
        disagreement origins
        add_finding keycloak-origin "Keycloak public URLs must share one origin: $REPLY"
    fi
    if [ "${#realms[@]}" -gt 1 ]; then
        disagreement realms
        add_finding keycloak-realm "Keycloak realm must be identical everywhere: $REPLY"
    fi
}

# When Apache is part of the generated topology, Keycloak's public URL must
# reach a VirtualHost that proxies to the Keycloak service. The service name is
# deliberately matched by its conventional keycloak token so applications may
# retain their own prefixes and network aliases.
check_keycloak_vhost() {
    local value host index name route_host matched=0
    local U_SCHEME U_NETLOC U_HOST U_PORT U_PORT_VALID U_PATH
    local -a names

    [ -n "${ENV[KEYCLOAK_PUBLIC_URL]+set}" ] || return 0
    trim "${ENV[KEYCLOAK_PUBLIC_URL]}"
    value="$REPLY"
    [ -n "$value" ] && [[ "$value" != *"$PLACEHOLDER"* ]] || return 0
    [ "${#ROUTE_TARGETS[@]}" -gt 0 ] || return 0

    url_split "$value"
    host="$U_HOST"
    [ -n "$host" ] && [ "$U_PORT_VALID" = 1 ] || return 0

    for index in "${!ROUTE_TARGETS[@]}"; do
        mapfile -t names <<< "${ROUTE_NAMES[$index]}"
        for name in "${names[@]}"; do
            [ -n "$name" ] || continue
            host_name "$name"
            [ "$REPLY" = "$host" ] || continue
            matched=1
            url_split "${ROUTE_TARGETS[$index]}"
            if [[ "$U_HOST" == *keycloak* ]]; then return 0; fi
            resolve_service "$U_HOST" || true
            if [[ "$REPLY" == *keycloak* ]]; then return 0; fi
        done
    done

    py_repr "$host"
    if [ "$matched" = 1 ]; then
        add_finding keycloak-vhost "KEYCLOAK_PUBLIC_URL host $REPLY must proxy to the Keycloak Compose service"
    else
        add_finding keycloak-vhost "KEYCLOAK_PUBLIC_URL host $REPLY has no Apache Keycloak VirtualHost"
    fi
}

strip_trailing_slashes() {
    REPLY="$1"
    while [[ "$REPLY" == */ ]]; do REPLY="${REPLY%/}"; done
}

check_canonical_urls() {
    local key value auth_key auth_value value_shown
    for key in "${ENV_KEYS[@]}"; do
        [[ "$key" == *NEXTAUTH_URL ]] || continue
        value="${ENV[$key]}"
        auth_key="${key%NEXTAUTH_URL}AUTH_URL"
        auth_value="${ENV[$auth_key]:-}"
        if [ -z "$value" ] || [ -z "$auth_value" ] \
                || [[ "$value$auth_value" == *"$PLACEHOLDER"* ]]; then
            continue
        fi
        strip_trailing_slashes "$value"
        value_shown="$REPLY"
        strip_trailing_slashes "$auth_value"
        if [ "$value_shown" != "$REPLY" ]; then
            py_repr "$value"
            value_shown="$REPLY"
            py_repr "$auth_value"
            add_finding canonical-url "$auth_key=$REPLY and $key=$value_shown must be identical"
        fi
    done
}

run_checks() {
    local mode="$1"
    apply_application_ports
    [ "$mode" != "production" ] || check_placeholders
    check_internal_urls
    check_django_hosts
    # Test Compose files set DB_HOST to the generated test database directly.
    [ "$mode" != "production" ] || check_databases
    check_keycloak
    check_keycloak_vhost
    check_canonical_urls
}

main() {
    local mode="" env_file="" compose_file="" apache_dir="" option value
    local name profile errors=0 warnings=0 level index
    while [ "$#" -gt 0 ]; do
        option="$1"
        case "$option" in
            -h|--help) print_help; exit 0 ;;
            --*=*) value="${option#*=}"; option="${option%%=*}"; shift ;;
            --mode|--env-file|--compose-file|--apache-config-dir|--service)
                [ "$#" -ge 2 ] || usage_error "argument $option: expected one argument"
                value="$2"
                shift 2
                ;;
            *) usage_error "unrecognized arguments: $option" ;;
        esac
        case "$option" in
            --mode) mode="$value" ;;
            --env-file) env_file="$value" ;;
            --compose-file) compose_file="$value" ;;
            --apache-config-dir) apache_dir="$value" ;;
            --service)
                name="${value%%=*}"
                profile="${value#*=}"
                if [[ "$value" != *=* ]] || [ -z "$name" ] || [ -z "$profile" ]; then
                    usage_error "invalid --service '$value'; expected NAME=PROFILE"
                fi
                [ -n "${PROFILES[$name]+set}" ] || PROFILE_ORDER+=("$name")
                PROFILES["$name"]="$profile"
                ;;
            *) usage_error "unrecognized arguments: $option" ;;
        esac
    done
    case "$mode" in
        production|test) ;;
        "") usage_error "the following arguments are required: --mode" ;;
        *) usage_error "argument --mode: invalid choice: '$mode' (choose from 'production', 'test')" ;;
    esac
    [ -n "$env_file" ] || usage_error "the following arguments are required: --env-file"
    [ -n "$compose_file" ] || usage_error "the following arguments are required: --compose-file"

    for value in "$env_file" "$compose_file"; do
        if [ ! -f "$value" ] || [ ! -r "$value" ]; then
            echo "Error: cannot read $value" >&2
            exit 2
        fi
    done

    COMPOSE_FILE="$compose_file"
    read_environment "$env_file"
    read_compose_services "$compose_file"
    if [ -n "$apache_dir" ] && [ -d "$apache_dir" ]; then
        read_apache_routes "$apache_dir"
    fi

    run_checks "$mode"
    for index in "${!FINDING_RULES[@]}"; do
        if [ "${FINDING_CONFIRMED[$index]}" = 1 ] && [ "$mode" = "production" ]; then
            level="ERROR"
            errors=$((errors + 1))
        else
            level="WARNING"
            warnings=$((warnings + 1))
        fi
        printf '%s [%s] %s\n' "$level" "${FINDING_RULES[$index]}" "${FINDING_MESSAGES[$index]}"
    done
    echo "Deployment configuration check ($mode): $errors error(s), $warnings warning(s)."
    [ "$errors" -eq 0 ]
}

main "$@"
