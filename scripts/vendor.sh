#!/usr/bin/env bash
# Copies nf-core/rnaseq and nf-core/differentialabundance into pipelines/nf-core/ and generates the
# config that an including pipeline has to provide itself (included pipelines do not bring their config).
#
# Override the sources with local checkouts while developing:
#   RNASEQ_SRC=~/projects/rnaseq-composable-input DIFFAB_SRC=~/projects/differentialabundance-composable scripts/vendor.sh
set -euo pipefail
cd "$(dirname "$0")/.."

RNASEQ_REPO=https://github.com/nf-core/rnaseq
RNASEQ_REF=21d0b581f          # feat/composable-input (nf-core/rnaseq#1966)
DIFFAB_REPO=https://github.com/nf-core/differentialabundance
DIFFAB_REF=0e60a7b1           # dev

vendor() { # name repo ref src-override
    local name=$1 repo=$2 ref=$3 src=${4:-}
    rm -rf "pipelines/nf-core/$name"; mkdir -p "pipelines/nf-core/$name"
    if [ -n "$src" ]; then
        rsync -a --exclude .git --exclude .nextflow --exclude work --exclude .nf-test "$src"/ "pipelines/nf-core/$name"/
    else
        local tmp; tmp=$(mktemp -d)
        git clone --quiet "$repo" "$tmp"; git -C "$tmp" checkout --quiet "$ref"
        rsync -a --exclude .git "$tmp"/ "pipelines/nf-core/$name"/; rm -rf "$tmp"
    fi
}

vendor rnaseq "$RNASEQ_REPO" "$RNASEQ_REF" "${RNASEQ_SRC:-}"
vendor differentialabundance "$DIFFAB_REPO" "$DIFFAB_REF" "${DIFFAB_SRC:-}"

mkdir -p conf/generated
RNA=pipelines/nf-core/rnaseq
DA=pipelines/nf-core/differentialabundance

# Config params of both pipelines. rnaseq's process config and differentialabundance's workflow read these
# as top-level params, so they have to exist here even though rnaseq also receives them as params.rnaseq.
sed -n '/^params {/,/^}/p' "$RNA/nextflow.config" | grep -v "outdir" > conf/generated/rnaseq_params.config
sed -n '/^params {/,/^}/p' "$DA/nextflow.config" > conf/generated/diffabundance_params.config

# Process config of the two pipelines. includeConfig resolves relative to the including file, so the
# generated file lists paths relative to conf/generated/. The selectors match the include aliases.
{
    grep -E "^includeConfig '(\./)?(conf/modules/|subworkflows/)" "$RNA/nextflow.config" \
        | sed -E "s#includeConfig '(\./)?#includeConfig '../../$RNA/#"
    echo "includeConfig '../../$DA/conf/modules.config'"
} > conf/generated/process.config
echo "vendored rnaseq and differentialabundance; generated conf/generated/*"
