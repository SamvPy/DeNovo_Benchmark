include { SPECTRALIS as SPECTRALIS_EA        } from "./refinement"
include { SPECTRALIS as SPECTRALIS_RESCORING } from "./refinement"
include { SPECTRALIS_PARSER                  } from "./refinement"
include { INSTANOVO_PLUS                     } from "./refinement"
include { INSTANOVO_PLUS_PARSER              } from "./refinement"
include { matchFiles                         } from "../utils/utility_fn"

workflow spectralis_rescoring_workflow {
    take:
        engine
        mgf_files
        result_files

    main:
        mgf_result_map = matchFiles(mgf_files, result_files)
        mgf_annotated = SPECTRALIS_PARSER(mgf_result_map, engine)


        mgf_annotated_fixed = mgf_annotated.map { item ->
            def (engine, mgf_output) = item
            def mgf_list = (mgf_output instanceof List) ? mgf_output : [mgf_output]
            tuple(engine, mgf_list)
        }

        rescoring_input = mgf_annotated_fixed.flatMap { engine, mgf_list -> 
            mgf_list.collect { mgf -> tuple( mgf, engine, 'rescoring') }
        }

        SPECTRALIS_RESCORING(rescoring_input, 'rescoring')
          
    emit:
        spectralis_processed = SPECTRALIS_RESCORING.out
}

workflow spectralis_ea_workflow {
    take:
        engine
        mgf_files
        result_files

    main:
        mgf_result_map = matchFiles(mgf_files, result_files)
        mgf_annotated = SPECTRALIS_PARSER(mgf_result_map, engine)

        rescoring_input = mgf_annotated.flatMap { engine, mgf_list -> 
            mgf_list.collect { mgf -> tuple( mgf, engine, 'ea') }
        }
        SPECTRALIS_EA(rescoring_input, 'ea')
             
    
    emit:
        spectralis_processed = SPECTRALIS_EA.out
}

workflow spectralis_workflow {
    take:
        engine
        mgf_files
        result_files

    main:
        if (params.spectralis_mode == 'ea') {
            spectralis_ea_workflow(engine, mgf_files, result_files)
        }
        else if (params.spectralis_mode == 'rescoring') {
            spectralis_rescoring_workflow(engine, mgf_files, result_files)
        }

    emit:
        spectralis_processed = (params.spectralis_mode == 'ea')            \
            ? spectralis_ea_workflow.out.spectralis_processed               \
            : spectralis_rescoring_workflow.out.spectralis_processed
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
        spectralis_rescoring = params.run_spectralis                 \
            ? spectralis_workflow.out.spectralis_processed                    \
            : null

        instanovoplus_output = params.run_instanovoplus              \
            ? instanovoplus_workflow.out.instanovoplus_output                     \
            : null
}