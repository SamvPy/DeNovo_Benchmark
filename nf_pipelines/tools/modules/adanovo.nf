process ADANOVO {
    conda "${params.conda_env_dir}/adanovo_env"
    maxForks params.maxforks_tool
    tag "Running AdaNovo on ${mgf_file.baseName}..."

    // Store the casanovo result in provided folder
    storeDir "${params.outdir_raw}/adanovo"

    input:
        path mgf_file
        path serializer
        path config_adanovo

    output:
        path "${mgf_file.baseName}.adanovo.mztab"
        path serializer

    script:
        """
        python -m \\
            AdaNovo.adanovo \\
                --mode=denovo \\
                --model=${params.model_path_adanovo} \\
                --peak_path=${mgf_file} \\
                --config=${config_adanovo} \\
                --output=${mgf_file.baseName}.adanovo.mztab
        """
}