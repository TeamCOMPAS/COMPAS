#!/usr/bin/env bash

set -euo pipefail

HERE="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"

export LD_LIBRARY_PATH="$HERE/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

exec "$HERE/bin/COMPAS" "$@"
