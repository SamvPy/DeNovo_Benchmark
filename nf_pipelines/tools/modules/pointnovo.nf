process POINTNOVO {
    conda "${params.conda_env_dir}/pointnovo_env"
    maxForks params.maxforks_tool
    tag "Running PointNovo on ${mgf_file.baseName}..."

    // Store the pointnovo result in provided folder
    storeDir "${params.outdir_raw}/pointnovo"

    input:
        path mgf_file
        path serializer
    output:
        path "${mgf_file.baseName}.pointnovo.tsv"
        path serializer

    script:
        """
        echo 'not implemented'
        """
}