#!/usr/bin/env bash
set -euo pipefail
IN="$1"; OUT="$2"
source "/mnt/d/CADD paper exercise/gnina/activate.sh" >/dev/null
obabel "$IN" -O "$OUT" 2>/dev/null
