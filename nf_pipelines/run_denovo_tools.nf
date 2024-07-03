process CASANOVO {
    conda '/home/sam/miniconda3/envs/casanovo_env'
    maxForks 1
    tag "Running CasaNovo search on ${mgf_file.baseName}..."

    // Store the casanovo result in provided folder
    publishDir "${params.denovo_results_dir}/casanovo", mode: "copy", saveAs: { filename ->
        "${mgf_file.baseName}.mztab"}, pattern: "*.mztab"
    
    // Store the created log file
    publishDir "${params.denovo_results_dir}/casanovo", mode: "copy", saveAs: { filename ->
        "${mgf_file.baseName}.log"}, pattern: "*.log"
    

    input:
        path mgf_file
        path config_casanovo
        path serializer

    output:
        path "${mgf_file.baseName}.mztab"
        path "${mgf_file.baseName}.log"
        path serializer

    script:
        """
        casanovo sequence -o ${mgf_file.baseName}.mztab --config=$config_casanovo --model=${params.model_path_casanovo} $mgf_file
        """
}


// Implement an error handling step (default GPU, retry on CPU) but implement in the script part (either in python or groovy)
process INSTANOVO {
    conda '/home/sam/miniconda3/envs/instanovo'
    maxForks 1
    tag "Running InstaNovo search on ${mgf_file.baseName}..."

    // Store the instanovo result in provided folder
    publishDir "${params.denovo_results_dir}/instanovo", mode: "copy", saveAs: { filename ->
        "${mgf_file.baseName}.csv"}, pattern: "*.csv"

    input:
        path mgf_file
        path config_instanovo
        path serializer

    output:
        path "${mgf_file.baseName}.csv"
        path serializer

    script:
        """
        python -m instanovo.utils.convert_to_ipc $mgf_file ${mgf_file.baseName}.ipc --source_type=mgf
        python -m instanovo.transformer.predict ${mgf_file.baseName}.ipc ${params.model_path_instanovo} --denovo --output_path=${mgf_file.baseName}.csv
        """
}


// process GRAPHNOVO {
//     conda '/home/sam/miniconda3/envs/genova'
// }


process CONTRANOVO {
    conda '/home/sam/miniconda3/envs/ContraNovo'
    maxForks 7
    tag "Running ContraNovo search on ${mgf_file.baseName}..."

    // Store the ContraNovo results
    publishDir "${params.denovo_results_dir}/contranovo", mode: "copy", saveAs: { filename ->
        "${mgf_file.baseName}.mztab"}, pattern: "*.mztab"
    publishDir "${params.denovo_results_dir}/contranovo", mode: "copy", saveAs: { filename ->
        "${mgf_file.baseName}.log"}, pattern: "*.log"

    input:
        path mgf_file
        path config_contranovo
        path serializer

    output:
        path "${mgf_file.baseName}.mztab"
        path "${mgf_file.baseName}.log"
        path serializer

    script:
        """
        python -m ContraNovo.ContraNovo --mode=denovo --config=$config_contranovo --peak_path=$mgf_file --model=$params.model_path_contranovo --output=${mgf_file.baseName}
        """

}


process NOVOB {
    conda '/home/sam/miniconda3/envs/NovoB'
    maxForks 1
    tag "Running NovoB search on ${mgf_file.baseName}..."

    // Store the novob results
    publishDir "${params.denovo_results_dir}/novob", mode: "copy", saveAs: { filename ->
        "${mgf_file.baseName}.tsv"}, pattern: "*.txt"

    input:
        path mgf_file
        path serializer

    output:
        path "result.txt"
        path serializer

    script:
        """
        python -m NovoB.Prediction -m ${params.model_path_novob} -i $mgf_file -b 16 -g
        """
}


process PEPNET {
    conda '/home/sam/miniconda3/envs/PepNet_env'
    maxForks 1
    tag "Running PepNet search on ${mgf_file.baseName}..."

    // Store the PepNet results
    publishDir "${params.denovo_results_dir}/pepnet", mode: "copy", saveAs: { filename ->
        "${mgf_file.baseName}.tsv"}, pattern: "*.tsv"
    

    input:
        path mgf_file
        path serializer
    
    output:
        path "${mgf_file.baseName}.tsv"
        path serializer
    
    script:
        """
        python -m PepNet.denovo --input $mgf_file --model ${params.model_path_pepnet} --output ${mgf_file.baseName}.tsv
        """
}


// process SPECTRALIS {
//     conda '/home/sam/miniconda3/envs/spectralis_env'
// }


// process INSTANOVO_PLUS {
//     conda '/home/sam/miniconda3/envs/instanovo'
// }

workflow {

    // if (params.run_graphnovo) {

    // }
    
    serializer = Channel.fromPath(params.serializer)

    if (params.run_contranovo) {
        mgf_files = Channel.fromPath(params.mgf_files)
        config_contranovo = Channel.fromPath(params.config_contranovo)

        if (params.serialize) {
            (contranovo_result, serializer) = CONTRANOVO(mgf_files, config_contranovo.first(), serializer.first())
        }
        else {
            (contranovo_result, _) = CONTRANOVO(mgf_files, config_contranovo.first(), serializer.first())
        }
        
    }


    if (params.run_novob) {
        mgf_files = Channel.fromPath(params.mgf_files)

        if (params.serialize) {
            (novob_result, serializer) = NOVOB(mgf_files, serializer.first())
        }
        else {
            (novob_result, _) = NOVOB(mgf_files, serializer.first())
        }
    }


    if (params.run_pepnet) {
        mgf_files = Channel.fromPath(params.mgf_files)

        if (params.serialize) {
            (pepnet_result, serializer) = PEPNET(mgf_files, serializer.first())
        }
        else {
            (pepnet_result, _) = PEPNET(mgf_files, serializer.first())
        }
    }


    if (params.run_casanovo) {
        mgf_files = Channel.fromPath(params.mgf_files)
        config_casanovo = Channel.fromPath(params.config_casanovo)

        if (params.serialize) {
            (casanovo_result, _, serializer) = CASANOVO(mgf_files, config_casanovo.first(), serializer.first())
        }
        else {
            (casanovo_result, _, _) = CASANOVO(mgf_files, config_casanovo.first(), serializer.first())
        }
    }

    
    if (params.run_instanovo) {
        mgf_files = Channel.fromPath(params.mgf_files)
        config_instanovo = Channel.fromPath(params.config_instanovo)

        if (params.serialize) {
            (instanovo_result, serializer) = INSTANOVO(mgf_files, config_instanovo.first(), serializer.first())
        }
        else{
            (instanovo_result, _) = INSTANOVO(mgf_files, config_instanovo.first(), serializer.first())
        }
    }

}