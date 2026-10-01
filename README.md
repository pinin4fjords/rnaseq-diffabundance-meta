# rnaseq + differentialabundance meta-pipeline

nf-core/rnaseq and nf-core/differentialabundance composed into one Nextflow pipeline with pipeline
composition ([nextflow-io/nextflow#7213](https://github.com/nextflow-io/nextflow/pull/7213)). Both
pipelines run in one session, one DAG and one work directory, so `-resume` covers the whole analysis and
differentialabundance starts as soon as rnaseq's merged gene matrices exist.

```
samples -> NFCORE_RNASEQ -> merged counts / lengths / GTF -> DIFFERENTIALABUNDANCE -> report
```

## Design choices

Goal: a meta-pipeline that is just `main.nf`, a config shell and a params file on top of the vendored pipelines, with
nothing lost. That needs the component pipelines to follow the
[pipeline composition ADR](https://github.com/nextflow-io/nextflow/blob/master/adr/20260608-pipeline-composition.md), so the
component branches ([rnaseq#1966](https://github.com/nf-core/rnaseq/pull/1966),
[differentialabundance#758](https://github.com/nf-core/differentialabundance/pull/758)) change the pipelines instead of
working around them here.

- **Included as pipelines.** The typed `params {}` block is the interface: one record per pipeline (`rnaseq.*`,
  `diffab.*`), no name collisions, outputs as channels. rnaseq's matrices feed differentialabundance directly: one DAG,
  one `-resume`.
- **Config is the meta's job, so the pipelines ship it as loadable files.** An included pipeline's `nextflow.config` is
  not loaded. Each pipeline keeps its params defaults and process config in their own files, included from its
  `nextflow.config` (standalone runs are unchanged) and from the meta, which adds only the manifest, plugin, executor,
  containers, resources and env.
- **Params are read in one place.** In a composed run the global `params` is the meta's, so code outside the entry
  workflow that reads it gets the wrong values. The pipelines pass their params record down, use `moduleDir` instead of
  `projectDir`, and use module templates instead of `bin/`, which is not on an included pipeline's process path.
- **Tool arguments are process inputs, not config closures.** `ext.args = { params.extra_star_align_args }` reads the
  meta's `params` when included, so `--rnaseq.extra_star_align_args` never arrives and the meta would need duplicate
  top-level params. As the ADR and Ben Sherman's answer on [#7213](https://github.com/nextflow-io/nextflow/pull/7213) say,
  defaults belong in the process and values are passed in. The argument policy moves into functions next to the call,
  the value is an optional field of the module's input record, and `task.ext.args` is still appended last as the
  override hook. The cost is that vendored modules differ from nf-core/modules, so this lives on the component
  branches.
- **Selectors work both ways.** Included processes get the include alias as a prefix, so selectors accept one.
- **Outputs.** A pipeline returns channels and the meta chooses what to publish. rnaseq returns `gene_quant` and `gtf`
  without publishing them.
- **Temporary.** A Nextflow build with #7213, and a patched nf-schema (`validateParameters` blocks on `Channel` params).

## Status

Experimental. It depends on unmerged work: the rnaseq changes in
[nf-core/rnaseq#1966](https://github.com/nf-core/rnaseq/pull/1966) (on top of the output-records work in #1945),
the differentialabundance change in [nf-core/differentialabundance#758](https://github.com/nf-core/differentialabundance/pull/758),
and a patched nf-schema. `pipelines.json` records the branch and commit of each pipeline.

## Requirements

- A Nextflow build that includes #7213 (not in a release yet), as `$NEXTFLOW`
- Docker
- A patched nf-schema (`2.7.2-channel.4`, loaded through `NXF_PLUGINS_TEST_REPOSITORY` by `scripts/run.sh`),
  because `validateParameters` blocks on `Channel` params in the released plugin

## Run

```bash
scripts/vendor.sh              # (re)copies the pipelines listed in pipelines.json
NEXTFLOW=/path/to/nextflow scripts/run.sh -stub-run   # wiring check
NEXTFLOW=/path/to/nextflow scripts/run.sh             # small test run on nf-core test data
```

The pipelines are committed under `pipelines/nf-core/`, so the project runs as cloned. `pipelines.json` has the same
shape as an nf-core `modules.json`: the repository, branch and commit each pipeline was copied from. To move a pipeline,
change its entry (or set `RNASEQ_REF` / `DIFFERENTIALABUNDANCE_REF` to a branch or commit) and rerun
`scripts/vendor.sh`, which writes the resolved commit back. `RNASEQ_SRC` and `DIFFAB_SRC` vendor local checkouts
instead, without touching `pipelines.json`. The copies leave out each pipeline's tests, CI files and docs (other than
`docs/images/`).

## Layout

- `main.nf` includes both pipelines, wires rnaseq's outputs into differentialabundance and declares the
  outputs that are published. The meta-pipeline declares its own `params` and `output` blocks; included
  pipelines contribute neither.
- differentialabundance is included like rnaseq (`params as DiffabParams`, `workflow as NFCORE_DIFFERENTIALABUNDANCE`):
  its options are the nested record `params.diffab`, and rnaseq's merged matrices, the annotation and the static
  sample metadata and contrasts are given as its files.
- `nextflow.config` is the configuration shell. An included pipeline contributes only its scripts, so the
  manifest, resources, container settings and the pipelines' own config are provided here: it includes rnaseq's
  `conf/params.config` (its config params only) and `conf/modules.config`, differentialabundance's
  `conf/modules.config`, sets how published files are written, and sets the process environment and shell that the
  pipelines' `nextflow.config` files set.
- `assets/` holds the test samplesheet, sample metadata and contrasts.

## Result

A real run on the nf-core test data (pseudo-alignment only, `scripts/run.sh` with the default `params.json`)
takes about two minutes and publishes `out/rnaseq/multiqc`, `out/rnaseq/quant` (merged gene counts, lengths
and TPM) and `out/differentialabundance/report`.

## Notes

- Each pipeline's options are set under its own name: `--rnaseq.input` is the rnaseq samplesheet (Nextflow loads it into one record per row, as for rnaseq's own `--input`), `--diffab.input` the sample metadata and `--diffab.contrasts` the contrasts, and every other option of a pipeline is available the same way, for example `--rnaseq.extra_star_align_args` or `--diffab.deseq2_alpha`.
- differentialabundance builds the paramset of the run from `params.diffab` and validates it against its own
  `nextflow_schema.json`. 
  The `rnaseq` profile of differentialabundance is a set of param values, which `params.json` gives as `diffab.*`.
- Use a fresh work directory (no `-resume` from another one): rnaseq only publishes files under the
  current work directory.
- rnaseq returns `gene_quant`, the merged gene matrices of its primary quantifier (alignment-based unless
  `rnaseq.skip_alignment` is set), and `gtf`, the reference annotation. The default `params.json` runs the
  pseudo-alignment path; the aligned path (STAR and Salmon, `rnaseq.skip_alignment` false, QC on) has also been run
  on the test data.
