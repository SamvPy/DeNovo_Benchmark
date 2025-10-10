process RESCORING_MS2RESCORE {
    conda "${params.conda_env_dir}/denovo_analysis_env"
    maxForks 1

    script:
        """
        python ${params.ms2rescore_script} \\
            --config=${params.ms2rescore_config}
        """
}

workflow {
    RESCORING_MS2RESCORE()
}