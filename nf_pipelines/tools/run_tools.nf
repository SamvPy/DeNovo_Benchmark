include { CASANOVO    } from "./modules/casanovo"
include { INSTANOVO   } from "./modules/instanovo"
include { CONTRANOVO  } from "./modules/contranovo"
include { PEPNET      } from "./modules/pepnet"
include { NOVOB       } from "./modules/novob"
include { PIPRIMENOVO } from "./modules/piprimenovo"
include { PIHELIXNOVO } from "./modules/pihelixnovo"
include { ADANOVO     } from "./modules/adanovo"
include { SMSNET      } from "./modules/smsnet"       // Not implemented yet
include { DEEPNOVO    } from "./modules/deepnovo"     // Not implemented yet
include { POINTNOVO   } from "./modules/pointnovo"    // Not implemented yet
include { POWERNOVO   } from "./modules/powernovo"    // Not implemented yet


// This is coded in longform without function as dynamic process invocation didn't seem to work...
workflow RUN_TOOLS {

    take:
        mgf_files
        serializer

    main:
        def results_all = Channel.empty()

        // CONTRANOVO
        if (params.run_contranovo) {
            config_contranovo = Channel.fromPath(params.config_contranovo)

            if (params.serialize) {
                (contranovo_result, serializer) = CONTRANOVO(mgf_files, serializer, config_contranovo.first())
            }
            else {
                (contranovo_result, _) = CONTRANOVO(mgf_files, serializer, config_contranovo.first())
            }
            results_all = results_all.concat(contranovo_result)
            
        }

        // NOVOB
        if (params.run_novob) {
            
            if (params.serialize) {
                (novob_result, serializer) = NOVOB(mgf_files, serializer)
            }
            else {
                (novob_result, _) = NOVOB(mgf_files, serializer)
            }
            results_all = results_all.concat(novob_result)
        }

        // PEPNET
        if (params.run_pepnet) {
            
            if (params.serialize) {
                (pepnet_result, serializer) = PEPNET(mgf_files, serializer)
            }
            else {
                (pepnet_result, _) = PEPNET(mgf_files, serializer)
            }
            results_all = results_all.concat(pepnet_result)
        }

        // CASANOVO
        if (params.run_casanovo) {
            config_casanovo = Channel.fromPath(params.config_casanovo)

            if (params.serialize) {
                (casanovo_result, _, serializer) = CASANOVO(mgf_files, serializer, config_casanovo.first())
            }
            else {
                (casanovo_result, _, _) = CASANOVO(mgf_files, serializer, config_casanovo.first())
            }
            results_all = results_all.concat(casanovo_result)
        }

        // INSTANOVO
        if (params.run_instanovo) {
            config_instanovo = Channel.fromPath(params.config_instanovo)

            if (params.serialize) {
                (instanovo_result, serializer) = INSTANOVO(mgf_files, serializer, config_instanovo.first())
            }
            else {
                (instanovo_result, _) = INSTANOVO(mgf_files, serializer, config_instanovo.first())
            }
            results_all = results_all.concat(instanovo_result)
        }

        // SMSNET
        if (params.run_smsnet) {
            config_smsnet = Channel.fromPath(params.config_smsnet)

            if (params.serialize) {
                (smsnet_result, serializer) = SMSNET(mgf_files, serializer, config_smsnet.first())
            }
            else {
                (smsnet_result, _) = SMSNET(mgf_files, serializer, config_smsnet.first())
            }
            results_all = results_all.concat(smsnet_result)
        }

        // DEEPNOVO
        if (params.run_deepnovo) {

            if (params.serialize) {
                (deepnovo_result, serializer) = DEEPNOVO(mgf_files, serializer)
            }
            else {
                (deepnovo_result, _) = DEEPNOVO(mgf_files, serializer)
            }
            results_all = results_all.concat(deepnovo_result)
        }

        // PIHELIXNOVO
        if (params.run_pihelixnovo) {
            config_pihelixnovo = Channel.fromPath(params.config_pihelixnovo)

            if (params.serialize) {
                (pihelixnovo_result, serializer) = PIHELIXNOVO(mgf_files, serializer, config_pihelixnovo.first())
            }
            else {
                (pihelixnovo_result, _) = PIHELIXNOVO(mgf_files, serializer, config_pihelixnovo.first())
            }
            results_all = results_all.concat(pihelixnovo_result)
        }

        // PIPRIMENOVO
        if (params.run_piprimenovo) {
            config_piprimenovo = Channel.fromPath(params.config_piprimenovo)

            if (params.serialize) {
                (piprimenovo_result, serializer) = PIPRIMENOVO(mgf_files, serializer, config_piprimenovo.first())
            }
            else {
                (piprimenovo_result, _) = PIPRIMENOVO(mgf_files, serializer, config_piprimenovo.first())
            }
            results_all = results_all.concat(piprimenovo_result)
        }

        // AdaNovo
        if (params.run_adanovo) {
            config_adanovo = Channel.fromPath(params.config_adanovo)

            if (params.serialize) {
                (adanovo_result, serializer) = ADANOVO(mgf_files, serializer, config_adanovo.first())
            }
            else {
                (adanovo_result, _) = ADANOVO(mgf_files, serializer, config_adanovo.first())
            }
            results_all = results_all.concat(adanovo_result)
        }
    emit:
        result = results_all
}



workflow {

    mgf_ch = Channel.fromPath(params.mgf_files)
    serializer_ch = Channel.fromPath(params.serializer)

    results_all = RUN_TOOLS(mgf_ch, serializer_ch.first())
    results_all.collect().view()

}