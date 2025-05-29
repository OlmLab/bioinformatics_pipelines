/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THIS IS BROKEN AT THE MOMENT
----------------------------------------------------------------------------------------
*/

workflow help {
    // Read version from VERSION file
    def versionFile = file("${projectDir}/VERSION")
    def version = versionFile.exists() ? versionFile.text.trim() : "Unknown"

    log.info "\nOlm Lab nf-metgenomics-pipelines Pipeline Version: ${version}\n"
    log.info "Usage: nextflow run main.nf --entry <subworkflow or workflow> [options]\n"

    def mainFile = file("${projectDir}/main.nf").text
    def workflows = []
    def subworkflows = []

    // Iterate through each line in the main file
    mainFile.eachLine { line ->
        def trimmedLine = line.trim()
        if (trimmedLine.startsWith("workflow ")) {
            def parts = trimmedLine.split(/\s+/)
            if (parts.length > 1) {
                def workflowName = parts[1]
                if (workflowName.startsWith("wf_")) {
                    workflows << workflowName
                } else if (workflowName.startsWith("sw_")) {
                    subworkflows << workflowName
                }
            }
        }
    }

    log.info "Available Workflows:\n"
    workflows.each { workflowName ->
        log.info "  - ${workflowName}"
    }

    log.info "\nAvailable Subworkflows:\n"
    subworkflows.each { subworkflowName ->
        log.info "  - ${subworkflowName}"
    }
}
