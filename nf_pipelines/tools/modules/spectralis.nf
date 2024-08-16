process SPECTRALIS {
    conda "${params.conda_env_dir}/spectralis_env"
    maxForks "${params.maxforks_tool}"
    tag "Using Spectralis on ${mgf_file.baseName}"

    publishDir "${params.outdir_spectralis}", mode: "copy", saveAs: { filename ->
        "${mgf_file.baseName}_${params.spectralis_mode}.csv"}, pattern: "*.csv"
    
    input:
        path mgf_file

    output:
        path "${mgf_file.baseName}_${params.spectralis_mode}.csv"

    script:
        """
        spectralis \\
            --mode=${params.spectralis_mode} \\
            --input_path=$mgf_file \\
            --output_path=${mgf_file.baseName}_${params.spectralis_mode}.csv \\
            --config=${params.config_spectralis}
        """
}