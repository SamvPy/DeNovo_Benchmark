process CASANOVO {
    conda "${params.conda_env_dir}/casanovo_env"
    maxForks params.maxforks_tool
    tag "Running Casanovo on ${mgf_file.baseName}..."

    // Store the casanovo result in provided folder
    storeDir "${params.outdir_raw}/casanovo"

    input:
        path mgf_file
        path serializer
        path config_casanovo

    output:
        path "${mgf_file.baseName}.casanovo.mztab"
        path serializer

    script:
        """
        casanovo \\
            sequence \\
                -o ${mgf_file.baseName}.casanovo.mztab \\
                --config=$config_casanovo \\
                --model=${params.model_path_casanovo} \\
                $mgf_file
        """
}