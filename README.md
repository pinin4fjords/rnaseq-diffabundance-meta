# rnaseq + differentialabundance meta-pipeline

nf-core/rnaseq and nf-core/differentialabundance composed into one Nextflow pipeline with pipeline
composition ([nextflow-io/nextflow#7213](https://github.com/nextflow-io/nextflow/pull/7213)). Both
pipelines run in one session, one DAG and one work directory, so `-resume` covers the whole analysis and
differentialabundance starts as soon as rnaseq's merged gene matrices exist.

```
samples -> NFCORE_RNASEQ -> merged counts / lengths / GTF -> DIFFERENTIALABUNDANCE -> report
```

## Status

Experimental. It depends on unmerged work: the rnaseq changes in
[nf-core/rnaseq#1966](https://github.com/nf-core/rnaseq/pull/1966) (on top of the output-records work in #1945),
the differentialabundance change in [nf-core/differentialabundance#758](https://github.com/nf-core/differentialabundance/pull/758),
and a patched nf-schema. `scripts/vendor.sh` pins the commits from those branches.

## Requirements

- A Nextflow build that includes #7213 (not in a release yet), as `$NEXTFLOW`
- Docker
- A patched nf-schema (`2.7.2-channel.3`, loaded through `NXF_PLUGINS_TEST_REPOSITORY` by `scripts/run.sh`),
  because `validateParameters` blocks on `Channel` params in the released plugin

## Run

```bash
scripts/vendor.sh              # copies the two pipelines to pipelines/nf-core/ and generates conf/generated/
NEXTFLOW=/path/to/nextflow scripts/run.sh -stub-run   # wiring check
NEXTFLOW=/path/to/nextflow scripts/run.sh             # small test run on nf-core test data
```

`scripts/vendor.sh` takes `RNASEQ_SRC` and `DIFFAB_SRC` to vendor local checkouts instead of the pinned refs.

## Layout

- `main.nf` includes both pipelines, wires rnaseq's outputs into differentialabundance and declares the
  outputs that are published. The meta-pipeline declares its own `params` and `output` blocks; included
  pipelines contribute neither.
- `lib/diffab.nf` builds differentialabundance's paramset from the rnaseq outputs.
- `nextflow.config` is the configuration shell. An included pipeline contributes only its scripts, so the
  manifest, resources, container settings, config params and the `ext` settings in `conf/modules/` of both
  pipelines are supplied here (`conf/generated/` is produced by `scripts/vendor.sh`).
- `assets/` holds the test samplesheet, sample metadata and contrasts.

## Result

A real run on the nf-core test data (pseudo-alignment only, `scripts/run.sh` with the default `params.json`)
takes about two minutes and publishes `out/rnaseq/multiqc`, `out/rnaseq/quant` (merged gene counts, lengths
and TPM) and `out/differentialabundance/report`.

## Notes

- `--samples` is a samplesheet that Nextflow loads into one record per row, like rnaseq's own `--input`.
- `rnaseq.*` params are passed as a nested record. Params that exist only in rnaseq's `nextflow.config`
  (for example `umitools_bc_pattern`) are not part of that record and are set at the top level.
- differentialabundance validates its paramset against its own `nextflow_schema.json`, which refers to
  `assets/schema_*.json` relative to the project root; `scripts/vendor.sh` copies those files into `assets/`.
- Use a fresh work directory (no `-resume` from another one): rnaseq only publishes files under the
  current work directory.
- `--quantification pseudo|aligned` selects which merged matrices rnaseq hands on. Only `pseudo` has been run; the
  `aligned` path (STAR and Salmon) has not been exercised.
