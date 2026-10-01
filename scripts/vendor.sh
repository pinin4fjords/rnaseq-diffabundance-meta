#!/usr/bin/env bash
# Copies nf-core/rnaseq and nf-core/differentialabundance into pipelines/nf-core/ and generates the
# config that an including pipeline has to provide itself (included pipelines do not bring their config).
# The copies, the generated config and assets/schema_*.json are committed, so the project runs as cloned;
# rerun this script to update them. pipelines/nf-core/VENDORED.md records what was copied.
#
# Override the sources with local checkouts while developing:
#   RNASEQ_SRC=~/projects/rnaseq-composable-input DIFFAB_SRC=~/projects/differentialabundance-composable scripts/vendor.sh
set -euo pipefail
cd "$(dirname "$0")/.."

RNASEQ_REPO=https://github.com/nf-core/rnaseq
RNASEQ_REF=cd7194c14          # feat/composable-input (nf-core/rnaseq#1966)
DIFFAB_REPO=https://github.com/nf-core/differentialabundance
DIFFAB_REF=e2bb4e0f           # feat/pipeline-composition (nf-core/differentialabundance#758)

# Not needed to run the pipelines
EXCLUDES=(--exclude .git --exclude .github --exclude .devcontainer --exclude .gitignore --exclude .gitattributes
          --exclude .nextflow --exclude work --exclude .nf-test --exclude tests
          --include 'docs/' --include 'docs/images/***' --exclude 'docs/*'
          --exclude '*.nf.test' --exclude '*.nf.test.snap')

vendor() { # name repo ref src-override
    local name=$1 repo=$2 ref=$3 src=${4:-}
    rm -rf "pipelines/nf-core/$name"; mkdir -p "pipelines/nf-core/$name"
    if [ -n "$src" ]; then
        rsync -a "${EXCLUDES[@]}" "$src"/ "pipelines/nf-core/$name"/
        ref=$(git -C "$src" rev-parse --short=8 HEAD 2>/dev/null || echo local)
    else
        local tmp; tmp=$(mktemp -d)
        git clone --quiet "$repo" "$tmp"; git -C "$tmp" checkout --quiet "$ref"
        rsync -a "${EXCLUDES[@]}" "$tmp"/ "pipelines/nf-core/$name"/; rm -rf "$tmp"
    fi
    echo "- $name: $repo at $ref" >> pipelines/nf-core/VENDORED.md
}

mkdir -p pipelines/nf-core; : > pipelines/nf-core/VENDORED.md
vendor rnaseq "$RNASEQ_REPO" "$RNASEQ_REF" "${RNASEQ_SRC:-}"
vendor differentialabundance "$DIFFAB_REPO" "$DIFFAB_REF" "${DIFFAB_SRC:-}"

mkdir -p conf/generated

# nf-schema resolves the "schema" entries of nextflow_schema.json against the project root, which is
# this project when the pipeline is included, so the schemas of differentialabundance's input and
# contrasts params have to exist at the same relative paths here.
mkdir -p assets
cp pipelines/nf-core/differentialabundance/assets/schema_*.json assets/
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
