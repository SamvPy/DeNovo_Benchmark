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
                --output_path=${mgf_file.baseName}.instanovo.csv \\
                --knapsack_path=${params.knapsack_instanovo}
        """
}

process INSTANOVO_V1 {
    conda "${params.conda_env_dir}/instanovo_env"
    maxForks params.maxforks_tool
    tag "Running InstaNovo_v1 on ${mgf_file.baseName}..."

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
        instanovo transformer predict \\
            --data-path=$mgf_file \\
            --output-path=${mgf_file.baseName}.instanovo.csv \\
            --instanovo-model=${params.model_path_instanovo} \\
            --denovo
        deactivate
        """
}