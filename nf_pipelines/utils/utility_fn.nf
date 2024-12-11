#!/usr/bin/env nextflow

def matchFiles(peak_files, result_files) {

    mgf_map = peak_files
        .map { file -> 
        tuple (file.baseName, file)    
    }

    result_map = result_files
        .map { file -> 
        def baseName = file.baseName.split('\\.')[0]
        tuple (baseName, file)    
    }

    mgf_result_map = mgf_map
        .mix(result_map)
        .groupTuple()
        .filter{ it ->
            def (baseName, files) = it
            files[0] != null && files[1] != null
        }
        .map { it ->
            def (baseName, files) = it
            def peak_file = files.find { it.name.endsWith(".${params.raw_extension}")}
            def result_file = files.find { !it.name.endsWith(".${params.raw_extension}")}
            return [peak_file, result_file]
        }
    return mgf_result_map
}

def collectFiles(root_dir, engines) {
    def channels = Channel.empty()

    // Iterate over each engine and its corresponding directory
    engines.each { engine ->
        def engine_dir = "${root_dir}/${engine}"

        // Ensure the directory exists before proceeding
        if (!file(engine_dir).exists()) {
            log.warn "Directory ${engine_dir} does not exist!"
            return
        }

        channel_engine =  Channel.fromPath("$engine_dir/*")
        channels = channels.concat(channel_engine)
    }

    // Return the collected files
    return channels
}