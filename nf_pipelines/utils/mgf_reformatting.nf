process MGF_REFORMAT {
    conda "${params.conda_env_dir}/denovo_analysis_env"
    maxForks 3
    tag "Reformatting ${mgf_file.baseName}"

    publishDir "${params.mgf_reformatted_dir}", mode: 'copy', saveAs: { filename ->
        "${mgf_file.baseName}.mgf"}, pattern: "${mgf_file.baseName}_reformatted.mgf"

    publishDir "${params.log_out}", mode: 'copy', saveAs: { filename ->
        "${mgf_file.baseName}.log"}, pattern: "denovo_output_parsing.log"

    input:
        path mgf_file

    output:
        path "${mgf_file.baseName}_reformatted.mgf"
        path "denovo_output_parsing.log"

    script:
        """
        python -m \\
            denovo_utils.parsers.scripts.reformat_mgf \\
                -i $mgf_file \\
                -c ${params.max_charge} \\
                -p ${params.min_peak_length}
        """
}

workflow {
    mgf_ch = Channel.fromPath(params.mgf_files)
    (reformatted, log_files) = MGF_REFORMAT(mgf_ch)
}