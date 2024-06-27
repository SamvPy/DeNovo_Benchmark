#!/usr/bin/env nextflow

def matchFiles(mgf_files, result_files) {

    mgf_map = mgf_files
        .map { file -> 
        tuple (file.baseName, file)    
    }

    result_map = result_files
        .map { file -> 
        tuple (file.baseName, file)    
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
            def mgf_file = files.find { it.name.endsWith(".mgf")}
            def result_file = files.find { !it.name.endsWith(".mgf")}
            return [mgf_file, result_file]
        }
    return mgf_result_map
}