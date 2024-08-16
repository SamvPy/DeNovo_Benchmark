process PEPNET {
    conda "${params.conda_env_dir}/pepnet_env"
    maxForks params.maxforks_tool
    tag "Running PepNet on ${mgf_file.baseName}..."

    // Store the PepNet results
    storeDir "${params.outdir_raw}/pepnet"
    
    input:
        path mgf_file
        path serializer
    
    output:
        path "${mgf_file.baseName}.pepnet.tsv"
        path serializer
    
    script:
        """
        python -m \\
            PepNet.denovo \\
                --input $mgf_file \\
                --model ${params.model_path_pepnet} \\
                --output ${mgf_file.baseName}.pepnet.tsv
        """
}