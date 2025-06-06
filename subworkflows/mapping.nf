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
    bowtie2_to_sorted_bam;
    index_bowtie2;
} from "${projectDir}/modules/alignment"

include {
    get_sequences_from_sra;
} from "${projectDir}/modules/download"

include {add_prefix_to_fasta;
} from "${projectDir}/modules/binning"

include {parse_input_reads;
} from "${projectDir}/subworkflows/sylph"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow bowtie2_sw {
    // Parse the input
    parse_input_reads()
    parse_input_genome()

    // Determine input reads based on whether we have paired-end data
    inputs = parse_input_reads.out.has_paired_end.value ? 
        parse_input_reads.out.reads1.merge(parse_input_reads.out.reads2) : 
        parse_input_reads.out.reads1

    bowtie2_to_sorted_bam(
        parse_input_reads.out.sample_names,
        parse_input_genome.out.bt2_basename,
        inputs,
        parse_input_genome.out.bt2_index,
        false
    )

    emit:
        sample_names = bowtie2_to_sorted_bam.out.sample_name
        sorted_bam = bowtie2_to_sorted_bam.out.sorted_bam
        single_end = parse_input_reads.out.has_paired_end.value
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELPER WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow parse_input_genome {
    main:
        // Check that at least one required input is provided
        if (!params.genome_db && !params.fastas && !params.fastas_list && !params.bowtie2_index) {
            error """
            One of the following parameters must be provided:
            - genome_db: A FASTA file containing reference genomes
            - fastas: A list of FASTA file paths
            - fastas_list: A file listing paths to individual FASTA files
            - bowtie2_index: Pre-built Bowtie2 index files
            """
        }

        // Ensure only one input type is provided
        def provided_inputs = [params.genome_db, params.fastas, params.fastas_list, params.bowtie2_index].findAll { it != null }
        if (provided_inputs.size() > 1) {
            error "Only one input type should be provided (genome_db, fastas, fastas_list, or bowtie2_index)"
        }

        // Handle pre-built Bowtie2 index
        if (params.bowtie2_index) {
            // Check that the base path exists
            if (!file(params.bowtie2_index).exists()) {
                error "Bowtie2 index base path does not exist: ${params.bowtie2_index}"
            }
            bt2_index = Channel.fromPath("${params.bowtie2_index}.*")
        }

        // Handle single genome database file
        else if (params.genome_db) {
            genome_db = file(params.genome_db)
            if (!genome_db.exists()) {
                error "Genome database file does not exist: ${params.genome_db}"
            }
            index_bowtie2(genome_db, "genome_db")
            bt2_index = index_bowtie2.out.bowtie2_index_files
            bt2_basename = genome_db
        }

        // Handle list of FASTA file paths
        else if (params.fastas) {
            fasta_files = Channel.fromPath(params.fastas)

            // Add prefix if requested
            if (params.add_fasta_prefix) {
                add_prefix_to_fasta(fasta_files)
                fasta_files = add_prefix_to_fasta.out.prefixed_fasta
            }

            // Concatenate files into a single database
            concatenate_files(fasta_files, "genomes_db.fasta")
            genome_db = concatenate_files.out.concatenated_file
            index_bowtie2(genome_db, "genome_db")
            bt2_index = index_bowtie2.out.bowtie2_index_files
        }
        
        // Handle FASTA list file
        else {
            // Parse input table and validate
            if (!file(params.fastas_list).exists()) {
                error "FASTA list file does not exist: ${params.fastas_list}"
            }
            table = tableToDict(file("${params.fastas_list}"))
            if (!table.containsKey("fasta_files")) {
                error "FASTA list file must contain a 'fasta_files' column"
            }
            
            fasta_files = Channel.fromPath(table["fasta_files"].collect{t->file(t)})

            // Add prefix if requested
            if (params.add_fasta_prefix) {
                add_prefix_to_fasta(fasta_files)
                fasta_files = add_prefix_to_fasta.out.prefixed_fasta
            }

            // Concatenate files into a single database
            concatenate_files(fasta_files, "genomes_db.fasta")
            genome_db = concatenate_files.out.concatenated_file
            index_bowtie2(genome_db, "genome_db")
            bt2_index = index_bowtie2.out.bowtie2_index_files
            bt2_basename = genome_db
        }

    emit:
        bt2_index
        bt2_basename
}