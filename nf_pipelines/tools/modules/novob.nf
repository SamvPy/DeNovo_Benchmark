process NOVOB {
    conda "${params.conda_env_dir}/novob_env"
    maxForks params.maxforks_tool
    tag "Running NovoB on ${mgf_file.baseName}..."

    // Store the novob results
    storeDir "${params.outdir_raw}/novob"

    input:
        path mgf_file
        path serializer

    output:
        path "${mgf_file.baseName}.novob.tsv"
        path serializer

    script:
        """
        python -m \\
            NovoB.Prediction \\
                -m ${params.model_path_novob} \\
                -i $mgf_file \\
                -b 16 \\
                -g
        mv result.txt ${mgf_file.baseName}.novob.tsv
        """
}
