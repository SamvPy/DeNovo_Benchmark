process CONTRANOVO {
    conda "${params.conda_env_dir}/contranovo_env"
    maxForks params.maxforks_tool
    tag "Running ContraNovo on ${mgf_file.baseName}..."

    // Store the ContraNovo results
    storeDir "${params.outdir_raw}/contranovo"

    input:
        path mgf_file
        path serializer
        path config_contranovo

    output:
        path "${mgf_file.baseName}.contranovo.mztab"
        path serializer

    script:
        """
        python -m \\
            ContraNovo.ContraNovo \\
                --mode=denovo \\
                --config=$config_contranovo \\
                --peak_path=$mgf_file \\
                --model=$params.model_path_contranovo \\
                --output=${mgf_file.baseName}
        mv ${mgf_file.baseName}.mztab ${mgf_file.baseName}.contranovo.mztab
        """
}