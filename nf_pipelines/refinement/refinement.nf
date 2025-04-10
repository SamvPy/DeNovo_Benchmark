process INSTANOVO_PLUS_PARSER {
    conda "${params.conda_env_dir}/denovo_analysis_env"
    maxForks 1
    tag "Parsing InstaNovo output to InstaNovo+ input for ${mgf_file.baseName}..."

    input:
        tuple path(mgf_file), path(result_files)
        val engine

    output:
        path "${mgf_file.baseName}.feather"  // Parsed instanovo_plus input file
        path "${mgf_file.baseName}.json" // index-spectrum_id mapping file

    script:
        """
        python -m denovo_utils.parsers.scripts.instanovoplus_input \\
            -m $mgf_file \\
            -d ${engine}
        """

}

process INSTANOVO_PLUS {
    conda "${params.conda_env_dir}/instanovo_env"
    maxForks 1
    tag "Refining sequence predictions with InstaNovo+ for ${result_file.baseName}..."

    input:
        path result_file
        path mapping_file
        val engine

    output:
        path "${result_file.baseName}.csv"

    script:
        """
        python ${params.instanovo_plus_run_script} \\
            -i $result_file \\
            -m ${params.model_path_instanovo_diffusion} \\
            -c $mapping_file \\
            -o ${params.denovo_results_dir}/instanovoplus/${engine} \\
            -d ${params.gpu_device}
        """
}

process SPECTRALIS_PARSER {
    conda "${params.conda_env_dir}/denovo_analysis_env"
    maxForks 1
    tag "Creating annotated mgf-file(s) for Spectralis for ${mgf_file.baseName}"

    input:
        // path mgf_file
        tuple path(mgf_file), path(result_file)
        val engine

    output:
        tuple val(engine), path("*.spectralis.mgf")

    script:
        // denovo_engines = params.denovo_engines.join(' ')
        """
        python -m denovo_utils.parsers.scripts.spectralis_input \\
            -r $result_file \\
            -m $mgf_file \\
            -d ${engine}
        """
}

process SPECTRALIS {
    conda "${params.conda_env_dir}/spectralis_env"
    maxForks 1
    tag "Using Spectralis on ${mgf_file.baseName}"

    publishDir "${params.denovo_results_dir}/spectralis/${engine}", mode: "copy", saveAs: { filename ->
        "${mgf_file.baseName}_${mode}.csv"}, pattern: "${mgf_file.baseName}_${mode}.csv"
    
    input:
        tuple path(mgf_file), val(engine)
        val mode

    output:
        path "${mgf_file.baseName}_${mode}.csv"

    script:
        """
        spectralis \\
            --mode=${mode} \\
            --input_path=$mgf_file \\
            --output_path=${mgf_file.baseName}_${mode}.csv \\
            --config=${params.config_spectralis}
        """
}

