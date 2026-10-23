#!/usr/bin/env bash
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.."
exec python3 code/requirements/install_environments.py
