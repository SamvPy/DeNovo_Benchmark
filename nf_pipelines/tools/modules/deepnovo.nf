process DEEPNOVO {
    conda "${params.conda_env_dir}/deepnovo_env"
    maxForks params.maxforks_tool
    tag "Running DeepNovo on ${mgf_file.baseName}..."

    // Store the DeepNovo result in provided folder
    storeDir "${params.outdir_raw}/deepnovo"

    input:
        path mgf_file
        path serializer
    output:
        path "${mgf_file.baseName}.deepnovo.tsv"
        path serializer

    script:
        """
        python -m DeepNovo.deepnovo_main \\
            --input_file ${mgf_file} \\
            --knapsack_path ${params.knapsack_deepnovo} \\
            --output_filename ${mgf_file.baseName}.deepnovo.tsv \\
            --train_dir ${params.model_dir} \\
            --decode \\
            --beam_search
        """
}