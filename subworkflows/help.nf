/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THIS IS BROKEN AT THE MOMENT
----------------------------------------------------------------------------------------
*/

workflow help {
    // Read version from VERSION file
    def versionFile = file("${projectDir}/VERSION")
    def version = versionFile.exists() ? versionFile.text.trim() : "Unknown"

    def sections = [
        "global"         : "Global Parameters",
        "sw_sylph"       : "Sylph Subworkflow"
    ]

    def configFile = file("${projectDir}/nextflow.config").text
    def paramsBlock = (configFile =~ /params\s*\{([^}]*)\}/)[0][1]

    log.info "\nOlm Lab nf-metgenomics-pipelines Pipeline Version: ${version}\n"
    log.info "Usage: nextflow run main.nf --entry <subworkflow or workflow> [options]\n"

    sections.each { section, description ->
        log.info "[$section]: $description"
        
        def matches = []
        if (section == "global") {
            // Extract global parameters (those after "Global Parameters" comment)
            def globalPattern = ~/(?m)^\s*\/\/\s*Global Parameters\s*\n((?:[^\n]*\n)*?)(?=\/\/|$)/
            def globalMatch = (paramsBlock =~ globalPattern)
            if (globalMatch.find()) {
                def globalBlock = globalMatch.group(1)
                def paramPattern = ~/(?m)^\s*(\w+)\s*=\s*([^\/]+)(?:\s*\/\/\s*(.*))?$/
                matches = globalBlock.findAll(paramPattern)
            }
        } else {
            // Extract section parameters
            def sectionPattern = ~/(?m)^\s*\/\/\s*${section}\s*Parameters\s*\n((?:[^\n]*\n)*?)(?=\/\/|$)/
            def sectionMatch = (paramsBlock =~ sectionPattern)
            if (sectionMatch.find()) {
                def sectionBlock = sectionMatch.group(1)
                def paramPattern = ~/(?m)^\s*(\w+)\s*=\s*([^\/]+)(?:\s*\/\/\s*(.*))?$/
                matches = sectionBlock.findAll(paramPattern)
            }
        }

        matches.each { match ->
            def paramName = match[1]
            def defaultValue = match[2].trim()
            def comment = match[3] ? match[3].trim() : ""
            log.info "  --${section == "global" ? "" : section + "."}${paramName}: Default = ${defaultValue}${comment ? " | " + comment : ""}"
        }
        log.info ""
    }
}