/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MODULE IMPORTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { 
    tableToDict;
} from "${projectDir}/modules/files"

include {
    bowtie2_sw;
} from "${projectDir}/subworkflows/mapping"

include {
    get_sequences_from_sra;
} from "${projectDir}/modules/download"

include {
    profile_with_instrain;
    make_stb_file_instrain;
} from "${projectDir}/modules/strain"

include {find_genes_prodigal} from "${projectDir}/modules/genes"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow instrain_profile_wf {
    // Parse the reads input
    parse_input_reads()

    // Parse the genome input
    parse_input_genome()

    // Run inStrain
    profile_with_instrain(
        parse_input_reads.out.sorted_bam, 
        parse_input_genome.out.genome_db,
        parse_input_genome.out.stb,
        parse_input_genome.out.genes
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUB WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow parse_input_reads {
    main:
        // Check that either reads or bams are provided
        if (!params.input_file && !params.input_bams) {
            error "Please provide either reads or bam files using the input_file or input_bams parameter"
        }

        // Handle reads input
        if (params.input_file) {
            bowtie2_sw()
            sorted_bam = bowtie2_sw.out.sorted_bam
            sample_names = bowtie2_sw.out.sample_names
        }
        // Handle bam input
        else if (params.input_bams) {
            table = tableToDict(file("${params.input_bams}"))
            sorted_bam = Channel.fromPath(table["bam_files"].collect{t->file(t)})
            sample_names = Channel.from(table["sample_name"])
        }

    emit:
        sample_names
        sorted_bam
}

workflow parse_input_genome {
    main:
        // Handle genome database input
        if (params.genome_db && !params.stb) {
            error "When providing a genome database (genome_db), you must also provide an STB file (stb)"
        }
        if (!params.genome_db) {
            if (params.add_fasta_prefix) {
                add_prefix_to_fasta(params.fasta_file)
                fasta_file = add_prefix_to_fasta.out.prefixed_fasta
            }
            concatenate_files(params.fasta_file, "genomes_db.fasta")
            genome_db = concatenate_files.out.concatenated_file
        } else {
            genome_db = file(params.genome_db)
        }

        // Handle STB file input
        if (!params.stb) {
            make_stb_file_instrain(params.fasta_file, "genomes_db")
            stb = make_stb_file_instrain.out.stb_file
        } else {
            stb = file(params.stb)
        }

        // Handle genes input
        if (!params.genes_file) {
            find_genes_prodigal(genome_db)
            genes = find_genes_prodigal.out.genes_fna
        } else {
            genes = file(params.genes_file)
        }

    emit:
        genome_db
        stb
        genes
}