process PIHELIXNOVO {
    conda "${params.conda_env_dir}/pihelixnovo_env"
    maxForks params.maxforks_tool
    tag "Running pihelixnovo on ${mgf_file.baseName}..."

    // Store the pihelixnovo result in provided folder
    storeDir "${params.outdir_raw}/powernovo"

    input:
        path mgf_file
        path serializer
    output:
        path "${mgf_file.baseName}.pihelixnovo.tsv"
        path serializer

    script:
        """
        echo 'not implemented'
        """
}