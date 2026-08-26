#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REDITOOLS="${ROOT}/tools/REDItools2/src/cineca/reditools.py"

if [[ ! -f "${REDITOOLS}" ]]; then
    echo "ERROR: REDItools2 submodule not found: ${REDITOOLS}" >&2
    echo "Run: git submodule update --init --recursive" >&2
    exit 2
fi

export PYTHONNOUSERSITE=1
exec "${PYTHON:-python3}" "${REDITOOLS}" "$@"
