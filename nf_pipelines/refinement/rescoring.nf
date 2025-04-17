// Processes which take an MGF and a result file and rescores
// the results with the function of a specific de novo engine
include { matchFiles                         } from "../utils/utility_fn"

process RESCORING_CASANOVO {
    conda "${params.conda_env_dir}/casanovo_analysis_env"
    maxForks 1
    tag "Rescoring psmlist (${psmlist.baseName}) from raw (${peak.baseName}) with Casanovo."

    publishDir "${params.rescoring_out}/casanovo", mode: "copy", saveAs: { filename ->
        "${psmlist.baseName}_casanovo.${params.psmlist_filetype}"}, pattern: "results.parquet"
    

    input:
        tuple path(peak), path(psmlist)

    output:
        path "results.parquet"

    script:
        """
        python ${params.casanovo_rescore_script} \\
            --mzml=$peak \\
            --psmlist=$psmlist \\
            --path_model=${params.path_model_casanovo} \\
            --path_config=${params.path_config_casanovo} \\
            --modification_mapping=${params.modification_mapping} \\
            --psmlist_filetype=${params.psmlist_filetype}
        """
}

// process RESCORING_INSTANOVO {
//     input:

//     output:

//     script:
// }

// process RESCORING_CONTRANOVO {
//     input:

//     output:

//     script:
// }

// process RESCORING_NOVOB {
//     input:

//     output:

//     script:
// }

// process RESCORING_PEPNET {
//     input:

//     output:

//     script:
// }

// process RESCORING_DEEPNOVO {
//     input:

//     output:

//     script:
// }

// process RESCORING_POINTNOVO {
//     input:

//     output:

//     script:
// }

// process RESCORING_SMSNET {
//     input:

//     output:

//     script:
// }

workflow {
    peak_files = Channel.fromPath(params.peak_files)
    psm_lists = Channel.fromPath(params.psm_lists)

    peak_psmlist_map = matchFiles(peak_files, psm_lists)
    peak_psmlist_map.view()

    RESCORING_CASANOVO(peak_psmlist_map)
}