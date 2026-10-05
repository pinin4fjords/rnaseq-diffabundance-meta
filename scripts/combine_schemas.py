#!/usr/bin/env python3
"""Writes nextflow_schema.json for the meta-pipeline from the schemas of the vendored pipelines.

Seqera Platform builds its launch form from the nextflow_schema.json at the repository root. Each
pipeline's options are one nested object (`rnaseq`, `diffab`), matching the records of the `params`
block in main.nf. The form only renders sections at the top level of a schema, so each pipeline's
sections are merged into the properties of its object.

A schema parameter is kept only if the pipeline's typed `params` block declares it: the others are
config params, which are not part of the pipeline record. Options that main.nf sets itself are left out.
A parameter stays required only if the params block gives it no default: Platform reads defaults from
nextflow.config, never from main.nf, so the form would otherwise ask for values the pipeline already has.
Defaults under `${projectDir}` point into the pipeline's own directory: the launch form sends schema
defaults as param values, and `projectDir` is the meta-pipeline's root when the pipeline is included.
"""
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

PIPELINES = [
    # key in params, vendored pipeline, title, options that main.nf sets
    ('rnaseq', 'rnaseq', 'nf-core/rnaseq', set()),
    ('diffab', 'differentialabundance', 'nf-core/differentialabundance', {'matrix', 'feature_length_matrix', 'gtf'}),
]


def params_block(main_nf):
    """Maps each param that the params block declares to whether it has a default."""
    block = re.search(r'^params\s*\{\n(.*?)^\}', main_nf.read_text(), re.M | re.S)
    return {name: bool(default) for name, default in re.findall(r'^\s+(\w+)\s*:[^=\n]*(=)?', block.group(1), re.M)}


def pipeline_section(key, name, title, wired):
    pipeline_dir = ROOT / 'pipelines' / 'nf-core' / name
    schema = json.loads((pipeline_dir / 'nextflow_schema.json').read_text())
    declared = params_block(pipeline_dir / 'main.nf')

    properties, required = {}, []
    for ref in schema['allOf']:
        group = schema['$defs'][ref['$ref'].split('/')[-1]]
        for param, definition in group['properties'].items():
            if param not in declared or param in wired:
                continue
            if 'schema' in definition:
                definition = {**definition, 'schema': f"pipelines/nf-core/{name}/{definition['schema']}"}
            default = definition.get('default')
            if isinstance(default, str) and default.startswith('${projectDir}/'):
                definition = {**definition, 'default': default.replace('${projectDir}/', f'${{projectDir}}/pipelines/nf-core/{name}/', 1)}
            properties[param] = definition
        required += [p for p in group.get('required', []) if p in properties and not declared[p]]

    options = {
        'type': 'object',
        'title': f'{title} options',
        'description': schema['description'],
        'properties': properties,
    }
    if required:
        options['required'] = required
    return {
        'title': title,
        'type': 'object',
        'fa_icon': 'fas fa-project-diagram',
        'properties': {key: options},
    }


def main():
    defs = {key: pipeline_section(key, name, title, wired) for key, name, title, wired in PIPELINES}
    schema = {
        '$schema': 'https://json-schema.org/draft/2020-12/schema',
        '$id': 'https://raw.githubusercontent.com/pinin4fjords/rnaseq-diffabundance-meta/main/nextflow_schema.json',
        'title': 'rnaseq-diffabundance pipeline parameters',
        'description': 'nf-core/rnaseq composed with nf-core/differentialabundance',
        'type': 'object',
        '$defs': defs,
        'allOf': [{'$ref': f'#/$defs/{key}'} for key in defs],
    }
    (ROOT / 'nextflow_schema.json').write_text(json.dumps(schema, indent=4) + '\n')
    for key, section in defs.items():
        print(f"{key}: {len(section['properties'][key]['properties'])} options")


if __name__ == '__main__':
    main()
