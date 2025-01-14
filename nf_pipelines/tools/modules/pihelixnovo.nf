process PIHELIXNOVO {
    conda "${params.conda_env_dir}/pihelixnovo_env"
    maxForks params.maxforks_tool
    tag "Running pihelixnovo on ${mgf_file.baseName}..."

    // Store the pihelixnovo result in provided folder
    storeDir "${params.outdir_raw}/pihelixnovo"

    input:
        path mgf_file
        path serializer
        path config_pihelixnovo

    output:
        path "${mgf_file.baseName}.pihelixnovo.tsv"
        path serializer

    script:
        """
        python -m \\
            HelixNovo.main \\
                --mode=denovo \\
                --config=$config_pihelixnovo \\
                --gpu=0 \\
                --output=denovo.log \\
                --peak_path=$mgf_file \\
                --model=${params.model_path_pihelixnovo}
        mv denovo_denovo.txt ${mgf_file.baseName}.pihelixnovo.tsv
        """
}