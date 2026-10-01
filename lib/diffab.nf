//
// Builds the paramset that nf-core/differentialabundance takes from the merged rnaseq outputs.
// Untyped, because differentialabundance is, and because it reads its defaults from the top-level
// params that the meta-pipeline's config provides (pipelines/nf-core/differentialabundance/conf/params.config).
//

include { getDefaultConfigurations      } from '../pipelines/nf-core/differentialabundance/subworkflows/local/utils_nfcore_differentialabundance_pipeline/main'
include { validateConfigurations        } from '../pipelines/nf-core/differentialabundance/subworkflows/local/utils_nfcore_differentialabundance_pipeline/main'
include { addDifferentialRuntimeParams  } from '../pipelines/nf-core/differentialabundance/subworkflows/local/utils_nfcore_differentialabundance_pipeline/main'

def diffabParamset(matrix, lengths, gtf, input, contrasts) {
    // The meta-pipeline's own params are not differentialabundance params (and hold Path objects, which
    // cannot be serialised during validation)
    def defaults = getDefaultConfigurations().collect { paramset ->
        paramset.findAll { key, _value -> !(key in ['rnaseq', 'samples', 'sample_metadata']) }
    }
    def configured = defaults.collect { paramset ->
        // Paths are given as strings, as they are on the command line: the paramset is serialised to
        // JSON during validation, which cannot handle Path objects
        paramset + [ matrix: matrix.toUriString(), feature_length_matrix: lengths.toUriString(), gtf: gtf.toUriString(), input: input.toUriString(), contrasts: contrasts.toUriString() ]
    }
    def paramset = validateConfigurations(configured)
        .collect { it -> addDifferentialRuntimeParams(it) }
        .first()
    return [
        id: paramset.study_name,
        paramset_name: paramset.paramset_name,
        params: paramset.findAll { key, _value -> key != 'paramset_name' }
    ]
}
