include { matchFiles } from '../utils/utility_fn'

process INSTANOVO_PLUS_PARSER {
    conda "${params.conda_env_dir}/denovo_analysis_env"
    maxForks 1
    tag "Parsing InstaNovo output to InstaNovo+ input for ${mgf_file.baseName}..."

    input:
        tuple path(mgf_file), path(result_file)

    output:
        path "${mgf_file.baseName}.feather"  // Parsed instanovo_plus input file
        path "${mgf_file.baseName}.json" // index-spectrum_id mapping file

    script:
        """
        python -m denovo_utils.parsers.scripts.instanovoplus_input \\
            -r $result_file \\
            -m $mgf_file \\
            -d ${params.denovo_engine}
        """

}

process INSTANOVO_PLUS {
    conda "${params.conda_env_dir}/instanovo_env"
    maxForks 1
    tag "Refining sequence predictions with InstaNovo+ for ${result_file.baseName}..."

    input:
        path result_file
        path mapping_file

    output:
        path "${result_file.baseName}.csv"

    script:
        """
        python ${params.instanovo_plus_run_script} \\
            -i $result_file \\
            -m ${params.model_path_instanovo_diffusion} \\
            -c $mapping_file \\
            -o ${params.denovo_results_dir}/instanovoplus
            -d ${params.gpu_device}
        """
}

process SPECTRALIS_PARSER {
    conda "${params.conda_env_dir}/denovo_analysis_env"
    maxForks 1
    tag "Creating annotated mgf-file for Spectralis for ${mgf_file.baseName}"

    input:
        path mgf_file
        // tuple path(mgf_file), path(result_file)

    output:
        path "${mgf_file.baseName}_annotated.mgf"

    script:
        denovo_engines = params.denovo_engines.join(' ')
        """
        python -m denovo_utils.parsers.scripts.spectralis_input \\
            -r ${params.result_root_dir} \\
            -m $mgf_file \\
            -d ${denovo_engines}
        """
}

process SPECTRALIS {
    conda "${params.conda_env_dir}/spectralis_env"
    maxForks 1
    tag "Using Spectralis on ${mgf_file.baseName}"

    publishDir "${params.denovo_results_dir}/spectralis", mode: "copy", saveAs: { filename ->
        "${mgf_file.baseName}_${params.spectralis_mode}.csv"}, pattern: "*.csv"
    
    input:
        path mgf_file

    output:
        path "${mgf_file.baseName}_${params.spectralis_mode}.csv"

    script:
        """
        spectralis \\
            --mode=${params.spectralis_mode} \\
            --input_path=$mgf_file \\
            --output_path=${mgf_file.baseName}_${params.spectralis_mode}.csv \\
            --config=${params.config_spectralis}
        """
}

workflow {
    
    if (params.run_instanovoplus) {

        result_files = Channel.fromPath(params.result_files)
        mgf_files = Channel.fromPath(params.mgf_files)
        mgf_result_map = matchFiles(mgf_files, result_files)

        (result_file_parsed, mapping_file) = INSTANOVO_PLUS_PARSER(mgf_result_map)
        _ = INSTANOVO_PLUS(result_file_parsed, mapping_file)
    }

    if (params.run_spectralis) {

        // result_files = Channel.fromPath(params.result_files)
        mgf_files = Channel.fromPath(params.mgf_files)
        // mgf_result_map = matchFiles(mgf_files, result_files)

        // mgf_annotated = SPECTRALIS_PARSER(mgf_result_map)

        mgf_annotated = SPECTRALIS_PARSER(mgf_files)
        _ = SPECTRALIS(mgf_annotated)
    }
}