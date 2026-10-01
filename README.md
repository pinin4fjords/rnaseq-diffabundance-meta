# rnaseq + differentialabundance meta-pipeline

nf-core/rnaseq and nf-core/differentialabundance composed into one Nextflow pipeline with pipeline
composition ([nextflow-io/nextflow#7213](https://github.com/nextflow-io/nextflow/pull/7213)). Both
pipelines run in one session, one DAG and one work directory, so `-resume` covers the whole analysis and
differentialabundance starts as soon as rnaseq's merged gene matrices exist.

```
samples -> NFCORE_RNASEQ -> merged counts / lengths / GTF -> DIFFERENTIALABUNDANCE -> report
```

## Design choices

The aim is the smallest meta-pipeline that keeps everything the two pipelines can do: this repository adds a short
`main.nf`, a configuration shell and a params file on top of the vendored pipelines. That only works if the pipelines
follow the standards in the [pipeline composition ADR](https://github.com/nextflow-io/nextflow/blob/master/adr/20260608-pipeline-composition.md),
so the component branches ([rnaseq#1966](https://github.com/nf-core/rnaseq/pull/1966),
[differentialabundance#758](https://github.com/nf-core/differentialabundance/pull/758)) change the pipelines instead of
working around them here. The reasons, most important first:

**Pipelines are included as pipelines.** `include { params as RnaseqParams ; workflow as NFCORE_RNASEQ }` takes the typed
`params {}` block as the interface, so each pipeline's options are a record under its own name (`rnaseq.*`, `diffab.*`)
and cannot collide, and its published outputs are returned as channels. Because rnaseq's merged matrices are channels,
differentialabundance starts as soon as they exist, and the whole analysis is one DAG that `-resume` covers.

**An included pipeline brings its scripts, not its config.** The ADR is explicit that the meta-pipeline provides the
configuration, and that process config can be reused if it sits in its own file. Both pipelines therefore keep their
config params defaults in `conf/params.config` and their process config in its own file, included from their
`nextflow.config` at the same place as before (the resolved config of a standalone run is unchanged). The meta includes
those files instead of copying or generating them, and adds only what the ADR says it must own: the manifest, the
plugin, the executor, container and resource settings, and the process environment and shell.

**Pipeline params are read in one place.** In a composed run the global `params` belongs to the meta-pipeline, so any
code outside the entry workflow that reads `params.x` silently gets the wrong value (a composed run of the unmodified
pipeline reported dozens of undefined params). The component branches pass the pipeline's params record down explicitly,
use `moduleDir` instead of `projectDir`, and keep scripts in module templates instead of `bin/`, which is not on the
path of an included pipeline's processes (a composed run failed with `command not found` for a script in `bin/`).

**Tool arguments are passed to processes, not read by config.** nf-core pipelines usually write
`ext.args = { params.extra_star_align_args ?: '' }` in their config. That works when the pipeline runs directly, but when
it is included the closure sees the meta-pipeline's `params`, not the `rnaseq` record: an option set as
`--rnaseq.extra_star_align_args` never reaches it, the meta-pipeline would have to declare copies of such params at the
top level and keep them in sync, and one pipeline could not be called twice with different tool settings. The ADR's
guidance, and Ben Sherman's answer on the [#7213 review](https://github.com/nextflow-io/nextflow/pull/7213), is that
Nextflow does not know which config belongs to an included workflow, so defaults should be defined in the process
definition and values passed as process inputs. In rnaseq the argument policy therefore moves from the config closures
into plain functions next to the call, and the value travels as an optional field of the module's input record.
`task.ext.args` is still appended last, so a meta-pipeline or user can override any tool setting with an alias-qualified
selector (`withName: 'NFCORE_RNASEQ:.*:STAR_ALIGN' { ext.args = ... }`). Only the process that uses a value receives it,
so changing one setting invalidates only the tasks it affects. The cost is that the vendored modules' input records
differ from nf-core/modules, which is why this lives on the component branches and would be proposed upstream
separately. Params that only affect configuration (publish mode, institutional config, resource caps, container
options) stay in config, as the ADR recommends. This change is in progress on rnaseq#1966; until it lands, `main.nf`
declares the rnaseq params that its process config still reads.

**Selectors must work both ways.** An included pipeline's processes are prefixed with the include alias, so a selector
such as `NFCORE_RNASEQ:RNASEQ:...` never matches in a composed run. The pipelines' selectors accept a prefix, and a
composed run is checked for selectors that match nothing.

**Outputs.** A pipeline's outputs are channels returned to its caller, and the meta-pipeline decides what to publish.
rnaseq returns `gene_quant` (the merged gene matrices of its primary quantifier) and `gtf` alongside its other outputs;
those two are declared with `enabled false`, so they are returned but not published by rnaseq. differentialabundance is
a typed pipeline with a typed `output {}` block, and takes the matrix, feature lengths and annotation as dataflow values.

**Vendored and pinned.** Until there is a pipeline registry the pipelines are copied under `pipelines/nf-core/` and
`pipelines.json` records where each came from, in the shape of an nf-core `modules.json`.

**Temporary pieces.** Two things are only here until they are released upstream: a Nextflow build that contains #7213,
and a patched nf-schema, because `validateParameters` blocks forever when a pipeline declares a `Channel` param.

## Status

Experimental. It depends on unmerged work: the rnaseq changes in
[nf-core/rnaseq#1966](https://github.com/nf-core/rnaseq/pull/1966) (on top of the output-records work in #1945),
the differentialabundance change in [nf-core/differentialabundance#758](https://github.com/nf-core/differentialabundance/pull/758),
and a patched nf-schema. `pipelines.json` records the branch and commit of each pipeline.

## Requirements

- A Nextflow build that includes #7213 (not in a release yet), as `$NEXTFLOW`
- Docker
- A patched nf-schema (`2.7.2-channel.3`, loaded through `NXF_PLUGINS_TEST_REPOSITORY` by `scripts/run.sh`),
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
  `conf/params.config` (defaults of the config params that its process config reads as top-level params) and
  `conf/process.config`, differentialabundance's `conf/modules.config`, sets the three top-level params that
  differentialabundance's process config reads (`outdir`, `publish_dir_mode`, `shinyngs_deploy_to_shinyapps_io`), and
  sets the process environment and shell that the pipelines' `nextflow.config` files set.
- `assets/` holds the test samplesheet, sample metadata and contrasts.

## Result

A real run on the nf-core test data (pseudo-alignment only, `scripts/run.sh` with the default `params.json`)
takes about two minutes and publishes `out/rnaseq/multiqc`, `out/rnaseq/quant` (merged gene counts, lengths
and TPM) and `out/differentialabundance/report`.

## Notes

- Each pipeline's options are set under its own name: `--rnaseq.input` is the rnaseq samplesheet (Nextflow loads it into one record per row, as for rnaseq's own `--input`), `--diffab.input` the sample metadata and `--diffab.contrasts` the contrasts, and every other option of a pipeline is available the same way, for example `--rnaseq.extra_star_align_args` or `--diffab.deseq2_alpha`.
- `rnaseq.*` params are passed as a nested record. rnaseq's process config reads some params as top-level params:
  `aligner`, `gencode`, `with_umi`, `pseudo_aligner` and the other params declared at the top of `main.nf` are also
  params of the pipeline itself, so they are set at the top level and passed into the record. Params that only the
  process config reads (for example `umitools_bc_pattern` or `extra_star_align_args`) are not part of the record and
  are also set at the top level.
- differentialabundance builds the paramset of the run from `params.diffab` and validates it against its own
  `nextflow_schema.json`, so the params that schema requires (here `outdir`) are part of the record. It also refers
  to `assets/schema_*.json` relative to the project root; `scripts/vendor.sh` copies those files into `assets/`.
  The `rnaseq` profile of differentialabundance is a set of param values, which `params.json` gives as `diffab.*`.
- Use a fresh work directory (no `-resume` from another one): rnaseq only publishes files under the
  current work directory.
- rnaseq returns `gene_quant`, the merged gene matrices of its primary quantifier (alignment-based unless
  `rnaseq.skip_alignment` is set), and `gtf`, the reference annotation. The default `params.json` runs the
  pseudo-alignment path; the aligned path (STAR and Salmon, `rnaseq.skip_alignment` false, QC on) has also been run
  on the test data.
