def version = file('VERSION').exists() ? file('VERSION').text.trim() : "Unknown"

include { help } from './subworkflows/help.nf'

workflow {
    if (params.help) {
        help()
        exit 0
    }

    log.info "Running pipeline version ${version}..."
    // Other workflows or processes can be called here
}


