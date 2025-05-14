workflow help {
    // Read version from VERSION file
    def versionFile = file('VERSION')
    def version = versionFile.exists() ? versionFile.text.trim() : "Unknown"

    def sections = [
        "global"         : "Global Parameters",
        "analyze_data"   : "Analyze Data Subworkflow",
        "preprocess_data": "Preprocess Data Subworkflow",
        "summarize_results": "Summarize Results Subworkflow"
    ]

    def configFile = file('nextflow.config').text
    def pattern = ~/params(?:\.(\w+))?\.(\w+)\s*=\s*([^\/\/]+)\s*\/\/\s*(.*)/

    log.info "\nOlm Lab nf-metgenomics-piplines Pipeline Version: ${version}\n"
    log.info "Usage: nextflow run main.nf --entry <subworkflow or workflow> [options]\n"

    sections.each { section, description ->
        log.info "[$section]: $description"
        def matches = configFile.findAll(pattern).findAll { match ->
            match[1] == section || (section == "global" && !match[1])
        }

        matches.each { match ->
            def paramName = match[2]
            def defaultValue = match[3].trim()
            def comment = match[4]
            log.info "  --${section ? section + '.' : ''}${paramName}: Default = ${defaultValue} | ${comment}"
        }
        log.info ""
    }
}