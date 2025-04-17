include { refinement_workflow as  workflow_casanovo    } from './refinement_workflows'
include { refinement_workflow as  workflow_instanovo   } from './refinement_workflows'
include { refinement_workflow as  workflow_contranovo  } from './refinement_workflows'
include { refinement_workflow as  workflow_novob       } from './refinement_workflows'
include { refinement_workflow as  workflow_pepnet      } from './refinement_workflows'
include { refinement_workflow as  workflow_smsnet      } from './refinement_workflows'
include { refinement_workflow as  workflow_deepnovo    } from './refinement_workflows'
include { refinement_workflow as  workflow_pointnovo   } from './refinement_workflows'
include { refinement_workflow as  workflow_piprimenovo } from './refinement_workflows'
include { refinement_workflow as  workflow_pihelixnovo } from './refinement_workflows'
include { refinement_workflow as  workflow_adanovo     } from './refinement_workflows'
include { refinement_workflow as  workflow_parquet     } from './refinement_workflows'


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
    if ('piprimenovo' in params.denovo_engines) {
        piprimenovo_results = Channel.fromPath("${params.result_root_dir}/piprimenovo/*")
        workflow_piprimenovo('piprimenovo', mgf_files, piprimenovo_results)
    }
    if ('pihelixnovo' in params.denovo_engines) {
        pihelixnovo_results = Channel.fromPath("${params.result_root_dir}/pihelixnovo/*")
        workflow_pihelixnovo('pihelixnovo', mgf_files, pihelixnovo_results)
    }
    if ('adanovo' in params.denovo_engines) {
        adanovo_results = Channel.fromPath("${params.result_root_dir}/adanovo/*")
        workflow_adanovo('adanovo', mgf_files, adanovo_results)
    }
    if ('parquet' in params.denovo_engines) {
        parquet_results = Channel.fromPath("${params.result_root_dir}/parquet/*")
        workflow_parquet('parquet', mgf_files, parquet_results)
    }
}