#!/bin/bash
# Usage: bash code/reproduce.sh /path/to/python /path/to/new-output [--optimized]
# No automatic installation or network access. Stored certificates are untouched.
set -euo pipefail
if [ "$#" -lt 2 ]; then
  echo "Usage: $0 PYTHON NEW_OUTPUT_DIRECTORY [--optimized]" >&2
  exit 2
fi
PYTHON_BIN="$1"
OUTPUT_DIRECTORY="$2"
shift 2
exec "$PYTHON_BIN" -B "$(dirname "$0")/replay.py" --out "$OUTPUT_DIRECTORY" "$@"
