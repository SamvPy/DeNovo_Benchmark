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
        path "${mgf_file.baseName}.powernovo.csv"
        path serializer

    script:
        """
        mkdir results
        python -m \\
            powernovo.run \\
                ${mgf_file} \\
                --working_folder=${params.model_path_powernovo} \\
                --output_folder='results' \\
                --batch_size=64 \\
                --assembler=false \\
                --protein_inference=false \\
                --use_bert=true
        mv results/${mgf_file.baseName}_pw_score.csv ${mgf_file.baseName}.powernovo.csv
        """
}