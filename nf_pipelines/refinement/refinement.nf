process INSTANOVO_PLUS_PARSER {
    conda "${params.conda_env_dir}/denovo_analysis_env"
    maxForks 1
    tag "Parsing ${engine} output to InstaNovo+ input for ${mgf_file.baseName}..."

    input:
        tuple path(mgf_file), path(result_files)
        val engine

    output:
        path "${mgf_file.baseName}.ipc"  // Parsed instanovo_plus input file
        path "${mgf_file.baseName}.csv" // index-spectrum_id mapping file
        val engine

    script:
        """
        python -m denovo_utils.parsers.scripts.instanovoplus_input_v1 \\
            -m $mgf_file \\
            -r $result_files \\
            -d ${engine}

        # Ensure expected output files exist even if script produced nothing
        if [ ! -f "${mgf_file.baseName}.ipc" ]; then
            echo "No .ipc file produced — creating empty placeholder."
            touch "${mgf_file.baseName}.ipc"
        fi

        if [ ! -f "${mgf_file.baseName}.csv" ]; then
            echo "No .csv file produced — creating empty placeholder."
            touch "${mgf_file.baseName}.csv"
        fi
        """

}

process INSTANOVO_PLUS {
    conda "${params.conda_env_dir}/instanovo_env"
    maxForks 1
    tag "Refining sequence predictions with InstaNovo+ for ${input_file.baseName}..."

    publishDir "${params.denovo_results_dir}/instanovoplus/${engine}", mode: "copy", saveAs: { filename ->
        "${input_file.baseName}.csv"}, pattern: "${input_file.baseName}.csv"

    input:
        path input_file
        path init_predictions
        path config_instanovo
        val engine

    output:
        path "${input_file.baseName}.csv"

    script:
        """
        if [ ! -s "$init_predictions" ]; then
            echo "init_predictions file is empty — skipping InstaNovo+ prediction."
            touch "${input_file.baseName}.csv"
        else
        # Copy config files to the InstaNovo package directory, otherwise IN+ always takes preinstalled configs...
            cp -r ./$config_instanovo/ \
            ${params.conda_env_dir}/instanovo_env/lib/python3.12/site-packages/instanovo/

            instanovo diffusion predict \\
                --data-path "./$input_file" \\
                --output-path "./${input_file.baseName}.csv" \\
                instanovo_predictions_path="./$init_predictions" \\
                --config-path="./configs/inference_test" \\
                --config-name="instanovoplus.yaml"
        fi
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
        
        # Check if any .spectralis.mgf files were created
        if ! ls *.spectralis.mgf 1> /dev/null 2>&1; then
            echo "No output files generated — creating empty placeholder."
            touch ${mgf_file.baseName}_1.spectralis.mgf
        fi
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
        if [ ! -s "$mgf_file" ]; then
            echo "Input file is empty — skipping Spectralis and creating empty CSV."
            touch ${mgf_file.baseName}_${mode}.csv
        else
            spectralis \\
                --mode=${mode} \\
                --input_path=$mgf_file \\
                --output_path=${mgf_file.baseName}_${mode}.csv \\
                --config=${params.config_spectralis}
        fi
        """
}

process INSTANOVO_PLUS_PARSER_LEGACY {
    conda "${params.conda_env_dir}/denovo_analysis_env"
    maxForks 1
    tag "Parsing ${engine} output to InstaNovo+ input for ${mgf_file.baseName}..."

    input:
        tuple path(mgf_file), path(result_files)
        val engine

    output:
        path "${mgf_file.baseName}.ipc"  // Parsed instanovo_plus input file
        path "${mgf_file.baseName}.json" // index-spectrum_id mapping file

    script:
        """
        python -m denovo_utils.parsers.scripts.instanovoplus_input \\
            -m $mgf_file \\
            -d ${engine}
        """

}

process INSTANOVO_PLUS_LEGACY {
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
        mkdir -p ${params.denovo_results_dir}/instanovoplus/${engine}
        python ${params.instanovo_plus_run_script} \\
            -i $result_file \\
            -m ${params.model_path_instanovo_diffusion} \\
            -c $mapping_file \\
            -o ${params.denovo_results_dir}/instanovoplus/${engine} \\
            -d ${params.gpu_device}
        """
}

