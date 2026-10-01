#!/usr/bin/env bash
# Copies nf-core/rnaseq and nf-core/differentialabundance into pipelines/nf-core/. The copies and
# assets/schema_*.json are committed, so the project runs as cloned; rerun this script to update them.
#
# Override the sources with local checkouts while developing:
#   RNASEQ_SRC=~/projects/rnaseq-composable-input DIFFAB_SRC=~/projects/differentialabundance-composable scripts/vendor.sh
set -euo pipefail
cd "$(dirname "$0")/.."

# Sources and pinned commits are recorded in pipelines.json (the same shape as an nf-core modules.json).
# To move a pipeline to another commit, change its branch and git_sha there (or set <NAME>_REF, e.g.
# RNASEQ_REF=feat/composable-input, to resolve a ref and write the commit back) and rerun this script.

# Not needed to run the pipelines
EXCLUDES=(--exclude .git --exclude .github --exclude .devcontainer --exclude .gitignore --exclude .gitattributes
          --exclude .nextflow --exclude work --exclude .nf-test --exclude tests
          --include 'docs/' --include 'docs/images/***' --exclude 'docs/*')

pin() { # repo-url org name field
    python3 -c "import json,sys; print(json.load(open('pipelines.json'))['repos'][sys.argv[1]]['pipelines'][sys.argv[2]][sys.argv[3]][sys.argv[4]])" "$@"
}

record() { # repo-url org name branch sha
    python3 - "$@" <<'PY'
import json, sys
url, org, name, branch, sha = sys.argv[1:]
data = json.load(open('pipelines.json'))
entry = data['repos'][url]['pipelines'][org][name]
entry['branch'], entry['git_sha'] = branch, sha
json.dump(data, open('pipelines.json', 'w'), indent=2)
open('pipelines.json', 'a').write('\n')
PY
}

vendor() { # org name src-override
    local org=$1 name=$2 src=${3:-}
    local repo="https://github.com/$org/$name"
    local dest="pipelines/$org/$name"
    local ref_var; ref_var=$(echo "$name" | tr a-z A-Z)_REF
    local branch; branch=$(pin "$repo" "$org" "$name" branch)
    local sha; sha=$(pin "$repo" "$org" "$name" git_sha)
    rm -rf "$dest"; mkdir -p "$dest"
    if [ -n "$src" ]; then
        rsync -a "${EXCLUDES[@]}" "$src"/ "$dest"/
        echo "warning: $name copied from $src, pipelines.json not updated"
    else
        local tmp; tmp=$(mktemp -d)
        git clone --quiet "$repo" "$tmp"
        git -C "$tmp" checkout --quiet "${!ref_var:-$sha}"
        sha=$(git -C "$tmp" rev-parse HEAD)
        [ -n "${!ref_var:-}" ] && branch=${!ref_var}
        rsync -a "${EXCLUDES[@]}" "$tmp"/ "$dest"/; rm -rf "$tmp"
        record "$repo" "$org" "$name" "$branch" "$sha"
    fi
}

vendor nf-core rnaseq "${RNASEQ_SRC:-}"
vendor nf-core differentialabundance "${DIFFAB_SRC:-}"

# nf-schema resolves the "schema" entries of nextflow_schema.json against the project root, which is
# this project when the pipeline is included, so the schemas of differentialabundance's input and
# contrasts params have to exist at the same relative paths here.
mkdir -p assets
cp pipelines/nf-core/differentialabundance/assets/schema_*.json assets/
echo "vendored rnaseq and differentialabundance"
