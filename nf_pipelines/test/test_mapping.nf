include { matchFiles; collectFiles } from '../utils/utility_fn'

process TEST {
    maxForks 1

    input:
        tuple path(mgf_file), path(result_file)

    output:
        path mgf_file

    script:
        """
        echo $result_file
        """

}

workflow {


        result_files = collectFiles(params.root_dir, params.engines)
        mgf_result_map = matchFiles(mgf_files, result_files)

        mgf_files = Channel.fromPath(params.mgf_dir)
        
        
        TEST(mgf_result_map)
}