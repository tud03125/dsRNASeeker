#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GETA2I="${ROOT}/tools/SPRINT/utilities/getA2I.py"
SPRINT_ENV="${SPRINT_CONDA_ENV:-sprint_env}"

if [[ ! -f "${GETA2I}" ]]; then
    echo "ERROR: SPRINT getA2I.py not found: ${GETA2I}" >&2
    echo "Run: git submodule update --init --recursive" >&2
    exit 2
fi

exec conda run -n "${SPRINT_ENV}" \
    --no-capture-output python "${GETA2I}" "$@"
