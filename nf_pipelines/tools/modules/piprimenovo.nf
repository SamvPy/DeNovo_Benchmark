process PIPRIMENOVO {
    conda "${params.conda_env_dir}/piprimenovo_env"
    maxForks params.maxforks_tool
    tag "Running pi-PrimeNovo on ${mgf_file.baseName}..."

    // Store the piprimenovo result in provided folder
    storeDir "${params.outdir_raw}/piprimenovo"

    input:
        path mgf_file
        path serializer
        path config_primenovo

    output:
        path "${mgf_file.baseName}.piprimenovo.tsv"
        path serializer

    script:
        """
        python -m \\
            PrimeNovo.PrimeNovo \\
                --mode=denovo \\
                --config=$config_primenovo \\
                --peak_path=$mgf_file \\
                --model=$params.model_path_piprimenovo
        mv denovo.tsv ${mgf_file.baseName}.piprimenovo.tsv
        """
}