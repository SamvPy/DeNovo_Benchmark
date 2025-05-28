// Processes which take an MGF and a result file and rescores
// the results with the function of a specific de novo engine
include { matchFiles                         } from "../utils/utility_fn"

process RESCORING_MS2RESCORE {
    conda "${params.conda_env_dir}/casanovo_analysis_env"
    maxForks 1
    tag "Rescoring psmlist (${psmlist.baseName}) from raw (${peak.baseName}) with MS2Rescore."

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

workflow {
    peak_files = Channel.fromPath(params.peak_files)
    psm_lists = Channel.fromPath(params.psm_lists)

    peak_psmlist_map = matchFiles(peak_files, psm_lists)
    peak_psmlist_map.view()

    RESCORING_MS2RESCORE(peak_psmlist_map)
}