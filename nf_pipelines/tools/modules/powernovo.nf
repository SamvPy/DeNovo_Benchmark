process POWERNOVO {
    conda "${params.conda_env_dir}/powernovo_env"
    maxForks params.maxforks_tool
    tag "Running PowerNovo on ${mgf_file.baseName}..."

    // Store the Powernovo result in provided folder
    storeDir "${params.outdir_raw}/powernovo"

    input:
        path mgf_file
        path serializer
    output:
        path "${mgf_file.baseName}.powernovo.tsv"
        path serializer

    script:
        """
        echo 'not implemented'
        """
}