process SMSNET {
    conda "${params.conda_env_dir}/smsnet_env"
    maxForks params.maxforks_tool
    tag "Running SMSNet on ${mgf_file.baseName}..."

    // Store the smsnet result in provided folder
    storeDir "${params.outdir_raw}/smsnet"

    input:
        path mgf_file
        path serializer
        path config_smsnet

    output:
        path "_output/${mgf_file.baseName}"
        path serializer

    script:
        """
        python -m \\
            SMSNet.run \\
                --model_dir ${params.model_path_smsnet} \\
                --inference_input_file ${mgf_file} \\
                --config ${config_smsnet} \\
                --rescore
        """
}