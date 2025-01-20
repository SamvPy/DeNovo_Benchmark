include { SPECTRALIS as SPECTRALIS_EA        } from "./refinement"
include { SPECTRALIS as SPECTRALIS_RESCORING } from "./refinement"
include { SPECTRALIS_PARSER                  } from "./refinement"
include { INSTANOVO_PLUS                     } from "./refinement"
include { INSTANOVO_PLUS_PARSER              } from "./refinement"
include { matchFiles                         } from "../utils/utility_fn"

workflow spectralis_workflow {
    take:
        engine
        mgf_files
        result_files

    main:
        mgf_result_map = matchFiles(mgf_files, result_files)
        mgf_annotated = SPECTRALIS_PARSER(mgf_result_map, engine)
        SPECTRALIS_EA(mgf_annotated, engine, 'ea')
        // The EA also reports the original scores.
        // SPECTRALIS_RESCORING(mgf_annotated, engine, 'rescoring')
    
    emit:
        spectralis_ea = SPECTRALIS_EA.out
        // spectralis_rescoring = SPECTRALIS_RESCORING.out
}

workflow instanovoplus_workflow {
    take:
        engine
        mgf_files
        result_files

    main:
        mgf_result_map = matchFiles(mgf_files, result_files)
        (result_file_parsed, mapping_file) = INSTANOVO_PLUS_PARSER(mgf_result_map, engine)
        INSTANOVO_PLUS(result_file_parsed, mapping_file, engine)
    
    emit:
        instanovoplus_output = INSTANOVO_PLUS.out
}

workflow refinement_workflow {
    take:
        engine
        mgf_files
        result_files

    main:
        if (params.run_instanovoplus) {
            instanovoplus_workflow(engine, mgf_files, result_files)
        }
        if (params.run_spectralis) {
            spectralis_workflow(engine, mgf_files, result_files)
        }

    emit:
        // spectralis_ea = spectralis_workflow.out.spectralis_ea
        // spectralis_rescoring = spectralis_workflow.out.spectralis_rescoring
        instanovoplus_output = instanovoplus_workflow.out.instanovoplus_output
}