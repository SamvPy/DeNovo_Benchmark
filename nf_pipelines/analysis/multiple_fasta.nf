include { SAGE } from "${params.path_search_module}"
include { matchFiles } from '../utils/utility_fn'

process ADD_FASTA_TO_CONFIG {
    tag "Adding ${fasta.baseName} to search configuration."

    input:
        path fasta
        path search_config

    output:
        path fasta
        path "${search_config.baseName}_new.json"
        val "${params.search_dir}/${fasta.baseName}"

    script:
        """
        #! ${HOME}/miniconda3/envs/denovo_analysis_env/bin/python
        import json

        json_path = '${search_config}'
        json_path_out = '${search_config}'.split(".")[0] + '_new.json'
        fasta_path = '${fasta}'

        with open(json_path, 'r') as f:
            json_file = json.load(f)
            json_file['database']['fasta'] = fasta_path

        with open(json_path_out, 'w') as f:
            json.dump(json_file, f)
        print(f'write fasta field in {json_path_out}.')
        """
}

process UNPAIR {
    input: 
        tuple path(fasta), path(config), val(publish_dir), path(spectrum)
    output:
        path fasta
        path config
        val publish_dir
        path spectrum
 
    script:
        """
        echo
        """
}

process ADD_DENOVO_SEQUENCES {
    publishDir "${params.store_fasta}", mode: "copy", pattern: "${fasta.baseName}_denovo.fasta"

    input:
        tuple path(fasta), path(spectrum_file), path(denovo_file)
        path search_config

    output:
        path spectrum_file
        path "${search_config.baseName}_new.json"
        path "${fasta.baseName}_denovo.fasta"
        val "${params.search_dir}/${fasta.baseName}_denovo"

    script:
        """
        #! ${HOME}/miniconda3/envs/denovo_analysis_env/bin/python

        import pandas as pd
        import json
        import os
        from denovo_utils.parsers.io.fasta import FastaHandler

        # ADD DENOVO SEQUENCES TO THE FASTA
        path_denovo_file = '$denovo_file'

        path_fasta = '$fasta'
        fasta_basename = os.path.basename(path_fasta).split(".")[0]
        fasta_path = f'{fasta_basename}_denovo.fasta'

        path_search_config = '${search_config}'
        path_search_config_out = path_search_config.split('.')[0] + '_new.json'

        denovo_engine = '${params.denovo_engine}'


        denovo_psms = pd.read_csv(path_denovo_file)
        denovo_psms = denovo_psms[denovo_psms.source==denovo_engine]

        fasta = FastaHandler()
        fasta.read(path_fasta)
        fasta.add_denovo_sequences(denovo_psms)
        boolean_filter = [True]*len(fasta.dataframe)

        fasta.write(
            boolean_filter=boolean_filter,
            out_path=f'./{fasta_basename}_denovo.fasta'
        )
        print(f'Written new fasta file: {fasta_basename}_denovo.fasta')


        # Add new fasta to config file and dump it as a new json file.
        with open(path_search_config, 'r') as f:
            json_file = json.load(f)
            json_file['database']['fasta'] = fasta_path

        with open(path_search_config_out, 'w') as f:
            json.dump(json_file, f)
        print(f'write fasta field in {path_search_config_out}.')
        """
}

workflow {

    spectrum_files = Channel.fromPath(params.spectrum_files)
    sage_config = Channel.fromPath(params.sage_config).first()
    fasta = Channel.fromPath(params.fasta_dir)
    
    if (params.add_denovo) {
        denovo_files = Channel.fromPath(params.denovo_files)

        spectrum_denovo_map = matchFiles(spectrum_files, denovo_files)

        fasta_spectrum_denovo_ch = fasta.combine(spectrum_denovo_map)
        (spectrum_files_unpaired, config_unpaired, fasta_unpaired, publish_dir_unpaired) = ADD_DENOVO_SEQUENCES(fasta_spectrum_denovo_ch, sage_config)

    } else {
        
        // Add fasta one-by-one to the search config file
        (fasta, sage_config_completed, publish_dir) = ADD_FASTA_TO_CONFIG(fasta, sage_config)

        // Create paired channel
        fasta_config_publish_ch = fasta.merge(sage_config_completed, publish_dir)
        branched_channel = fasta_config_publish_ch.combine(spectrum_files)

        // Unpair to make it compatible with SAGE process
        (fasta_unpaired, config_unpaired, publish_dir_unpaired, spectrum_files_unpaired) = UNPAIR(branched_channel)
    }
    

    // Perform search
    (search_results, spectrum_files, configs, quant) = SAGE(spectrum_files_unpaired, config_unpaired, fasta_unpaired, publish_dir_unpaired)
}
