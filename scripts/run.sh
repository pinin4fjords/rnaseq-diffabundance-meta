#!/usr/bin/env bash
# Usage: scripts/run.sh [nextflow args]   e.g. scripts/run.sh -stub-run
# Needs a Nextflow build that includes pipeline composition (nextflow-io/nextflow#7213) as $NEXTFLOW.
set -euo pipefail
cd "$(dirname "$0")/.."
export NXF_PLUGINS_TEST_REPOSITORY=https://github.com/pinin4fjords/nf-schema/releases/download/2.7.2-channel.3/nf-schema-2.7.2-channel.3-meta.json
exec "${NEXTFLOW:-nextflow}" run main.nf -params-file params.json -output-dir out "$@"
