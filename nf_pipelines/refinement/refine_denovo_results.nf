include { refinement_workflow as  workflow_casanovo } from './refinement_workflows'
include { refinement_workflow as  workflow_instanovo } from './refinement_workflows'
include { refinement_workflow as  workflow_contranovo } from './refinement_workflows'
include { refinement_workflow as  workflow_novob } from './refinement_workflows'
include { refinement_workflow as  workflow_pepnet } from './refinement_workflows'
include { refinement_workflow as  workflow_smsnet } from './refinement_workflows'
include { refinement_workflow as  workflow_deepnovo } from './refinement_workflows'
include { refinement_workflow as  workflow_pointnovo } from './refinement_workflows'


// Main workflow: Loop over the list of engines
workflow {

    mgf_files = Channel.fromPath(params.mgf_files)

    if ('casanovo' in params.denovo_engines) {
        casanovo_results = Channel.fromPath("${params.result_root_dir}/casanovo/*")
        workflow_casanovo('casanovo', mgf_files, casanovo_results)
    }
    if ('instanovo' in params.denovo_engines) {
        instanovo_results = Channel.fromPath("${params.result_root_dir}/instanovo/*")
        workflow_instanovo('instanovo', mgf_files, instanovo_results)
    }
    if ('contranovo' in params.denovo_engines) {
        contranovo_results = Channel.fromPath("${params.result_root_dir}/contranovo/*")
        workflow_contranovo('contranovo', mgf_files, contranovo_results)
    }
    if ('novob' in params.denovo_engines) {
        novob_results = Channel.fromPath("${params.result_root_dir}/novob/*")
        workflow_novob('novob', mgf_files, novob_results)
    }
    if ('pepnet' in params.denovo_engines) {
        pepnet_results = Channel.fromPath("${params.result_root_dir}/pepnet/*")
        workflow_pepnet('pepnet', mgf_files, pepnet_results)
    }
    if ('deepnovo' in params.denovo_engines) {
        deepnovo_results = Channel.fromPath("${params.result_root_dir}/deepnovo/*")
        workflow_deepnovo('deepnovo', mgf_files, deepnovo_results)
    }
    if ('pointnovo' in params.denovo_engines) {
        pointnovo_results = Channel.fromPath("${params.result_root_dir}/pointnovo/*")
        workflow_pointnovo('pointnovo', mgf_files, pointnovo_results)
    }
    if ('smsnet' in params.denovo_engines) {
        smsnet_results = Channel.fromPath("${params.result_root_dir}/smsnet/*")
        workflow_smsnet('smsnet', mgf_files, smsnet_results)
    }
}