process INSTANOVO {
    conda "${params.conda_env_dir}/instanovo_env"
    maxForks params.maxforks_tool
    tag "Running InstaNovo on ${mgf_file.baseName}..."

    // Store the instanovo result in provided folder
    storeDir "${params.outdir_raw}/instanovo"

    input:
        path mgf_file
        path serializer
        path config_instanovo

    output:
        path "${mgf_file.baseName}.instanovo.csv"
        path serializer

    script:
        """
        python -m \\
            instanovo.utils.convert_to_ipc \\
                $mgf_file \\
                ${mgf_file.baseName}.ipc \\
                --source_type=mgf
        python -m \\
            instanovo.transformer.predict \\
                ${mgf_file.baseName}.ipc \\
                ${params.model_path_instanovo} \\
                --denovo \\
                --output_path=${mgf_file.baseName}.instanovo.csv
        """
}