// NOT FUNCTIONAL. This is a proposal of how I could have implemented the workflow
// with nice dynamic process invocations which make extensibility trivial

include { CASANOVO   } from "./modules/casanovo"
include { INSTANOVO  } from "./modules/instanovo"
include { CONTRANOVO } from "./modules/contranovo"
include { PEPNET     } from "./modules/pepnet"
include { NOVOB      } from "./modules/novob"

// Function to run a given tool
def run_tool(tool_name, mgf_ch, serializer, config_name = null) {

    // Some de novo tools require a configuration file to be passed
    if (config_name) {
        config_name = Channel.fromPath(params["config_${tool_name}"]).first()
    }

    def tool_process = this.&"${tool_name.toUpperCase()}"
    return tool_process(mgf_ch, serializer, config_name ? config_name: null)

}


// This is coded in longform without function as dynamic process invocation didn't seem to work...
workflow RUN_TOOLS {

    take:
        mgf_files
        serializer

    main:
        def results_all = Channel.empty()

        params.tools.each { tool_name, tool_run -> 

            if (tool_run) {
                def config_name = params.tool_configs.containsKey(tool_name) ? params.tool_configs[tool_name] : null
                
                if (params.serialize) {
                    (result, serializer) = run_tool(tool_name, serializer, config_name)
                } else {
                    (result, _) = run_tool(tool_name, serializer, config_name)
                }
                results_all = results_all.concat(result)
            }
        }

    emit:
        results_all
}



workflow {

    mgf_ch = Channel.fromPath(params.mgf_files)
    serializer_ch = Channel.fromPath(params.serializer)

    results_all = RUN_TOOLS(mgf_ch, serializer_ch.first())

}