//
// Builds the paramset that nf-core/differentialabundance takes from the merged rnaseq outputs.
// Untyped, because differentialabundance is, and because it reads its defaults from the top-level
// params that the meta-pipeline's config provides (conf/generated/diffabundance_params.config).
//

include { getDefaultConfigurations      } from '../pipelines/nf-core/differentialabundance/subworkflows/local/utils_nfcore_differentialabundance_pipeline/main'
include { validateConfigurations        } from '../pipelines/nf-core/differentialabundance/subworkflows/local/utils_nfcore_differentialabundance_pipeline/main'
include { addDifferentialRuntimeParams  } from '../pipelines/nf-core/differentialabundance/subworkflows/local/utils_nfcore_differentialabundance_pipeline/main'

def diffabParamset(matrix, lengths, gtf, input, contrasts) {
    // The nested rnaseq record and the samples channel are not differentialabundance params
    def defaults = getDefaultConfigurations().collect { paramset ->
        paramset.findAll { key, _value -> !(key in ['rnaseq', 'samples']) }
    }
    def configured = defaults.collect { paramset ->
        paramset + [ matrix: matrix, feature_length_matrix: lengths, gtf: gtf, input: input, contrasts: contrasts ]
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
