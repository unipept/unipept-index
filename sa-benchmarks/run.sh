#!/usr/bin/env bash
# Entry point for the benchmark suites. Everything below it is Python; this exists so the driver
# can be invoked from the repo root without knowing where the package lives or how it is spelled.
#
#   ./sa-benchmarks/run.sh defaults              the production-defaults sweep (regression gate)
#   ./sa-benchmarks/run.sh ram --dry-run         plan a sweep without touching the index
#   sudo ./sa-benchmarks/run.sh all              every suite, into one report.md
#
# See sa-benchmarks/README.md for the suites and the full-database runbook.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

command -v python3 >/dev/null || { echo "ERROR: python3 not found"; exit 1; }

# `tomllib` is stdlib from 3.11; the driver deliberately has no third-party dependencies so it can
# run on a benchmark server without a virtualenv.
python3 - <<'PY' || { echo "ERROR: python3 >= 3.11 required (tomllib)"; exit 1; }
import sys
sys.exit(0 if sys.version_info >= (3, 11) else 1)
PY

# PYTHONPATH rather than `cd "$HERE"`: the driver must stay in the caller's directory, or every
# relative path they pass (--out, --baseline) would silently resolve against sa-benchmarks/.
export PYTHONPATH="$HERE${PYTHONPATH:+:$PYTHONPATH}"
exec python3 -m bench "$@"
