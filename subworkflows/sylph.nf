/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MODULE IMPORTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { 
    tableToDict;
    concatenate_files;
} from "${projectDir}/modules/files"

include {
    estimate_abundance_sylph;
    estimate_abundance_sylph_SE;
} from "${projectDir}/modules/abundance"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow sylph_sw {
    // Parse the input
    parse_input()

    // Set up database
    sylph_db = file(params.sylph_db)

    // Run Sylph analysis on each sample
    if (parse_input.out.has_paired_end.value) {
        estimate_abundance_sylph(
            parse_input.out.reads1_sylph,
            parse_input.out.reads2_sylph,
            sylph_db
        )
    } else {
        estimate_abundance_sylph_SE(
            parse_input.out.reads1_sylph,
            sylph_db
        )
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELPER WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow parse_input {
    main:
        // Validate required parameters
        if (!params.input_file) {
            error "Required parameter 'input_file' is missing"
        }

        // Parse input table and validate columns
        def input_table = tableToDict(file("${params.input_file}"))
        def required_columns = ['reads1', 'sample_name']
        def missing_columns = required_columns.findAll { !input_table.containsKey(it) }
        
        if (missing_columns) {
            error """
            Missing required columns: ${missing_columns.join(', ')}
            Required columns: reads1, sample_name
            Optional columns: reads2 (for paired-end reads)
            Found columns: ${input_table.keySet().join(', ')}
            """
        }

        // Create channels for reads and sample names
        reads_1 = Channel.fromPath(input_table["reads1"].collect { file(it) })
        sample_names = Channel.fromList(input_table["sample_name"])
        
        // Handle paired-end reads if available
        has_paired_end = input_table.containsKey("reads2")
        if (!has_paired_end) {
            log.warn "⚠️  Single-end reads detected - paired-end reads recommended for better results"
        }
        
        reads_2 = has_paired_end ? Channel.fromPath(input_table["reads2"].collect { file(it) }) : Channel.empty()

        // Process reads for Sylph
        reads_1.multiMap { t -> reads_sylph: t }.set { reads1_all }
        
        if (has_paired_end) {
            reads_2.multiMap { t -> reads_sylph: t }.set { reads2_all }
            sylph_reads = [reads1_all.reads_sylph, reads2_all.reads_sylph]
        } else {
            sylph_reads = [reads1_all.reads_sylph, Channel.empty()]
        }

    emit:
        sample_names
        reads1_sylph = params.input_type == "sra" ? sylph_reads.reads1_sylph : sylph_reads[0]
        reads2_sylph = params.input_type == "sra" ? sylph_reads.reads2_sylph : sylph_reads[1]
        has_paired_end
}