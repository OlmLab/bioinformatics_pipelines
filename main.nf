#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INITIALIZATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def version = file("${projectDir}/VERSION").exists() ? file("${projectDir}/VERSION").text.trim() : "Unknown"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUB-WORKFLOW ACCESS POINT LOADING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { help      }               from './subworkflows/help.nf'
include { sylph_sw  }               from './subworkflows/sylph.nf'
include { bowtie2_sw  }             from './subworkflows/mapping.nf'

include { instrain_profile_wf  }    from './workflows/instrain.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEFAULT WORKFLOW WITH NO ENTRY SPECIFIED
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    if (params.help) {
        help()
        exit 0
    }

    log.info "Running pipeline version ${version}. Please specify an -entry to get started"
    // Other workflows or processes can be called here
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ENTRY WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow sw_sylph {
    log.info "Running Sylph with pipeline version ${version}"
    sylph_sw ()
}

workflow sw_bowtie2 {
    log.info "Running Bowtie2 with pipeline version ${version}"
    bowtie2_sw ()
}

workflow wf_instrain_profile {
    log.info "Running inStrain with pipeline version ${version}"
    instrain_profile_wf ()
}