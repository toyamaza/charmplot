#!/usr/bin/env bash
# bash tempalte taken from: https://betterdev.blog/minimal-safe-bash-script-template/

set -Eeuo pipefail
trap cleanup SIGINT SIGTERM ERR EXIT

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd -P)

usage() {
    cat <<EOF
Usage: $(basename "${BASH_SOURCE[0]}") [-h] [-v] [--dry-run] -o param_value arg1 arg2 [arg3...]

Merge TRExFitter output folders given by arg1, arg2, ... into a folder specified with the -o argument.

Available options:

-h, --help       Print this help and exit
-v, --verbose    Print script debug info
-o, --output     Output folder
--dry-run        Do not execute trex-fitter
EOF
    exit
}

cleanup() {
    trap - SIGINT SIGTERM ERR EXIT
    # script cleanup here
}

setup_colors() {
    if [[ -t 2 ]] && [[ -z "${NO_COLOR-}" ]] && [[ "${TERM-}" != "dumb" ]]; then
        NOFORMAT='\033[0m' RED='\033[0;31m' GREEN='\033[0;32m' ORANGE='\033[0;33m' BLUE='\033[0;34m' PURPLE='\033[0;35m' CYAN='\033[0;36m' YELLOW='\033[1;33m'
    else
        NOFORMAT='' RED='' GREEN='' ORANGE='' BLUE='' PURPLE='' CYAN='' YELLOW=''
    fi
}

msg() {
    echo >&2 -e "${1-}"
}

die() {
    local msg=$1
    local code=${2-1} # default exit status 1
    msg "$msg"
    exit "$code"
}

parse_params() {
    # default values of variables set from params
    output=''
    dry_run=0

    while :; do
        case "${1-}" in
        -h | --help) usage ;;
        -v | --verbose) set -x ;;
        --dry-run) dry_run=1 ;;
        -o | --output)
            output="${2-}"
            shift
            ;;
        -?*) die "Unknown option: $1" ;;
        *) break ;;
        esac
        shift
    done

    args=("$@")

    # check required params and arguments
    [[ -z "${output-}" ]] && die "Missing required parameter: 'output'"
    [[ ${#args[@]} -eq 1 ]] && die "Missing script arguments"

    return 0
}

parse_params "$@"
setup_colors

# script logic here

msg "Read parameters:${NOFORMAT}"
msg "- output: ${output}"
msg "- arguments:"
for line in ${args[@]}; do
    msg "  > ${line}"
done

# check that hadd exists

if ! type "hadd" >/dev/null; then
    die "Did not find the hadd command"
fi

# make output folder

msg "Creating output folder ${output}..."

mkdir -p ${output}

# merge

for file in $(ls ${args[0]}); do
    if [[ ${dry_run} -eq 0 ]]; then
        hadd ${output}/${file} $(for path in ${args[@]}; do echo -n "${path}/${file} "; done)
    else
        echo "${output}/${file}" $(for path in ${args[@]}; do echo -n "${path}/${file} "; done)
    fi
done
