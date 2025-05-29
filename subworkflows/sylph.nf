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
    estimate_abundance_sylph(
        parse_input.out.reads1_sylph,
        parse_input.out.reads2_sylph,
        sylph_db
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELPER WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow parse_input {
    main:
        // Check required parameters
        def required_params = ['input_file']
        def missing_params = required_params.findAll { params[it] == null }
        
        if (missing_params) {
            error """
            The following required parameters are missing:
            ${missing_params.join('\n')}
            
            Please provide all required parameters to run this workflow.
            """
        }

        if (params.input_type=="sra") {
            table = tableToDict(file("${params.input_file}"))
            get_sequences_from_sra(Channel.fromList(table["Run"]))
            sample_names = get_sequences_from_sra.out.sra_ids
            reads = get_sequences_from_sra.out.fastq_files
        }
        else {
            table = tableToDict(file("${params.input_file}"))
            
            // Check required columns exist in input table
            def required_columns = ['reads1', 'reads2', 'sample_name'] 
            def missing_columns = required_columns.findAll { !table.containsKey(it) }
            
            if (missing_columns) {
                error """
                The following required columns are missing from the input file:
                ${missing_columns.join('\n')}
                
                Input file must contain columns: reads1, reads2, and sample_name
                """
            }

            // Create channels for each sample's reads
            reads_1 = Channel.fromPath(table["reads1"].collect{t->file(t)})
            reads_2 = Channel.fromPath(table["reads2"].collect{t->file(t)})
            sample_names = Channel.fromList(table["sample_name"])
        }

        // Process reads for different analyses
        if (params.input_type=="sra") {
            reads.multiMap{t->
                reads_sylph: t
                reads_metaphlan: t
            }.set{reads_all}

            reads_all.reads_sylph.multiMap{t->
                reads1_sylph: t[0]
                reads2_sylph: t[1]
            }.set{sylph_reads}
        } else {
            // For local files, reads are already paired
            reads_1.multiMap{t->
                reads_sylph: t
                reads_metaphlan: t
            }.set{reads1_all}

            reads_2.multiMap{t->
                reads_sylph: t
                reads_metaphlan: t
            }.set{reads2_all}

            sylph_reads = [reads1_all.reads_sylph, reads2_all.reads_sylph]
        }

    emit:
        sample_names
        reads1_sylph = params.input_type=="sra" ? sylph_reads.reads1_sylph : sylph_reads[0]
        reads2_sylph = params.input_type=="sra" ? sylph_reads.reads2_sylph : sylph_reads[1]
        reads_metaphlan = params.input_type=="sra" ? reads_all.reads_metaphlan : reads1_all.reads_metaphlan
}