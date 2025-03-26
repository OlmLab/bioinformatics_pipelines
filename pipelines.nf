params.output_dir = "./output_is"
params.paired = false
params.host_indexed = false
params.binning_extension = "fa"
params.add_fasta_prefix = true
params.is_genome_db=null // default is null
params.is_stb_db=null // default is null

// ###### MAIN WORKFLOW ###### //
workflow {
    if (params.roadmap_id=="roadmap_1")
    {
        if (!params.host_genome)
            {
                error "Please provide a host genome for decontamination in the roadmap_1 workflow."
            }


            if (params.input_type=="sra")
            {
                table=tableToDict(file("${params.input_file}"))
                get_sequences_from_sra(Channel.fromList(table["Run"]))
                roadmap_1(get_sequences_from_sra.out.sra_ids, get_sequences_from_sra.out.fastq_files, file(params.host_genome))
            }
            if (params.input_type=="local")
            {
                table=tableToDict(file("${params.input_file}"))
                reads_1=Channel.fromPath(table["reads1"])
                reads_2=Channel.fromPath(table["reads2"])
                reads=reads_1.merge(reads_2)
                sample_name=Channel.fromList(table["sample_name"])
                roadmap_1(sample_name, reads, file(params.host_genome))
            }
            }
        
        else if (params.roadmap_id=="roadmap_2")
        {
            table=tableToDict(file("${params.input_reads}"))
            reads_1=Channel.fromPath(table["reads1"].collect{t->file(t)})
            reads_2=Channel.fromPath(table["reads2"].collect{t->file(t)})
            reads=reads_1.merge(reads_2)
            sample_names=Channel.from(table["sample_name"])
            
            if (params.input_fastas)
            {
            if (params.is_genome_db||params.is_stb_db)
            {
                error "If you provide a list of fasta files, you cannot provide a genome database or STB file."
            }
            else
            {
            fasta_file=tableToDict(file("${params.input_fastas}"))["fasta_files"].collect{t->file(t)}
            }
            }

            if (params.is_genome_db)
            {
                if (!params.is_stb_db)
                {
                    error "If you provide a genome database, you must also provide an STB file."
                }
                else{
                    fasta_file=file(params.is_genome_db)
                }
            }
            roadmap_2(sample_names, reads, fasta_file)
        }



        else if (params.roadmap_id=="roadmap_3")
        {
            if (params.input_type=="path")
            {
                genomes=files("${params.input_genomes}")

                roadmap_3(genomes)
            }
        }
        else if (params.roadmap_id=="roadmap_1_3_2")
        {
            if (!params.host_genome)
            {
                error "Please provide a host genome for decontamination in the roadmap_1 workflow."
            }
            if (params.input_type=="sra")
            {
                table=tableToDict(file("${params.input_file}"))
                get_sequences_from_sra(Channel.fromList(table["Run"]))
                sample_names=get_sequences_from_sra.out.sra_ids.collect()
                reads=get_sequences_from_sra.out.fastq_files.collect()
                roadmap_1_3_2(sample_names, reads, file(params.host_genome))
            }
            if (params.input_type=="local")
            {
                table=tableToDict(file("${params.input_file}"))
                reads_1=Channel.fromPath(table["reads1"])
                reads_2=Channel.fromPath(table["reads2"])
                reads=reads_1.merge(reads_2)
                sample_name=table["sample_name"]
                sample_name=Channel.fromList(sample_name)
                roadmap_1_3_2(sample_name, reads, file(params.host_genome))
            }
        }
    else if (params.roadmap_id=="roadmap_4")
    {
        if (!params.host_genome)
            {
                error "Please provide a host genome for decontamination in the roadmap_1 workflow."
            }


            if (params.input_type=="sra")
            {
                table=tableToDict(file("${params.input_file}"))
                get_sequences_from_sra(Channel.fromList(table["Run"]))
                roadmap_4(get_sequences_from_sra.out.sra_ids, get_sequences_from_sra.out.fastq_files, file(params.host_genome))
            }
            if (params.input_type=="local")
            {
                table=tableToDict(file("${params.input_file}"))
                reads_1=Channel.fromPath(table["reads1"])
                reads_2=Channel.fromPath(table["reads2"])
                reads=reads_1.merge(reads_2)
                sample_name=Channel.fromList(table["sample_name"])
                roadmap_4(sample_name, reads, file(params.host_genome))
            }
    }

    else if (params.roadmap_id=="roadmap_3_2")
    {
            if (!params.input_reads)
            {
                error "Please provide the reads information using the input_reads parameter."
            }
            if (!params.input_fastas)
            {
                error "Please provide the fasta files information using the input_fastas parameter."
            }
            table=tableToDict(file("${params.input_reads}"))
            reads_1=Channel.fromPath(table["reads1"].collect{t->file(t)})
            reads_2=Channel.fromPath(table["reads2"].collect{t->file(t)})
            reads=reads_1.merge(reads_2)
            sample_names=Channel.from(table["sample_name"])
            fasta_file=tableToDict(file("${params.input_fastas}"))["fasta_files"].collect{t->file(t)}
            roadmap_3_2(sample_names, reads, fasta_file)
    
    }

    else
        {
            error "Please provide a valid roadmap_id."
        }
        
        
    }
    

// ###### Roadmaps ###### //


workflow roadmap_1{
    // This roadmap takes the raw sequencing reads and performs:
    // 1- quality control
    // 2- assembly
    // 3- binning.

    take:
    sample_name
    reads
    host_genome

    main:
    quality_control(sample_name, reads, host_genome)
    assembly(quality_control.out.sample_name, quality_control.out.qc_reads)
    binning(assembly.out.sample_name, assembly.out.sorted_bams, assembly.out.contigs,assembly.out.reads)
    emit:
    metabat2_bins=binning.out.metabat2_bins
}


workflow roadmap_2 {
    // This roadmap takes a list of reads, their corresponding sample names, and a list of fasta files.
    // It performs strain-level profiling using InStrain.

    take:
    sample_names
    reads
    fasta_file
    main:
    if (!params.is_genome_db)
    {
        if (params.add_fasta_prefix)
        {
            add_prefix_to_fasta(fasta_file)
            fasta_file=add_prefix_to_fasta.out.prefixed_fasta
        }
        else
        {
            fasta_file=fasta_file
        }
        concatenate_files(fasta_file, "genomes_db.fasta")
        genome_db=concatenate_files.out.concatenated_file
    }
    else
    {
        genome_db=file(params.is_genome_db)
    }

    if (!params.is_stb_db)
    {

        make_stb_file_instrain(fasta_file, "genomes_db")
        stb=make_stb_file_instrain.out.stb_file
    }
    else
    {
        stb=file(params.is_stb_db)
    }

    if (!params.is_genes)
    {
        find_genes_prodigal(fasta_file)
        genes=find_genes_prodigal.out.genes_fna
    }
    else
    {
        genes=file(params.is_genes)
    }


    index_bowtie2(genome_db, "genomes_db")
    align_bowtie2(sample_names, index_bowtie2.out.reference_genome, reads, index_bowtie2.out.bowtie2_index_files)
    convert_sam_to_sorted_bam(align_bowtie2.out.bowtie2_sam, align_bowtie2.out.sample_name, align_bowtie2.out.paired)
    profile_with_instrain(convert_sam_to_sorted_bam.out.sorted_bam,genome_db, stb, genes)
    profile_with_instrain.out.instrain_profiles.collect().set{all_profiles}
    compare_instrain_profiles(all_profiles, stb)
    emit:
    instrain_profiles=all_profiles


}
    
workflow roadmap_3 {
    take:
    genomes
    
    main:
    write_genome_list(genomes)
    dereplicate_drep(genomes,write_genome_list.out.genomes_list)

    emit:
    dereplicated_genomes=dereplicate_drep.out.dereplicated_genomes

}

workflow roadmap_4 {
    take:
    sample_name
    reads
    host_genome
    main:
    quality_control(sample_name, reads, host_genome)
    emit:
    qc_reads=quality_control.out.qc_reads
    
}

workflow roadmap_1_3_2{
    /*
    This workflow gets the bins from individual samples using roadmap 1, dereplicates them using roadmap 3, and then performs strain-level profiling using roadmap 2.
    WARNING: This workflow is not fully tested and may not work as expected.
    */
    take:
    sample_name
    reads
    host_genome

    main:
    roadmap_1(sample_name, reads, host_genome)
    params.genomes_exctension=params.binning_extension
    dereplicated_genomes=roadmap_3(roadmap_1.out.metabat2_bins.collect())
    params.is_genome_db=null
    params.is_stb_db=null
    roadmap_2(sample_name, reads, dereplicated_genomes)


}

workflow roadmap_3_2 {
    // This roadmap is appropriate for ddereplicating a series of bins and then performing strain-level profiling using InStrain on the dereplicated bins.
    // An example would be comparative genomics for closely related bins
    take:
    sample_names
    reads
    fasta_file
    main:
    roadmap_3(fasta_file.collect())
    dereplicated_genomes=roadmap_3.out.dereplicated_genomes
    roadmap_2(sample_names, reads, dereplicated_genomes)

    emit:
    instrain_profiles=roadmap_2.out.instrain_profiles

}

// ###### WORKFLOWS ###### //

workflow quality_control {
    /*
    This workflow takes raw sequencing reads and performs quality control and decontamination.
    It generates cleaned reads and quality reports.
    */
    take:
    sample_name
    reads
    host_genome

    main:
    read_qc_fastp(sample_name, reads)
    index_bowtie2(host_genome,host_genome.baseName)
    align_bowtie2(read_qc_fastp.out.sample_name,index_bowtie2.out.reference_genome, read_qc_fastp.out.fastp_qcd_reads,index_bowtie2.out.bowtie2_index_files)  
    convert_sam_to_sorted_bam(align_bowtie2.out.bowtie2_sam,align_bowtie2.out.sample_name,align_bowtie2.out.paired)
    get_unmapped_reads(convert_sam_to_sorted_bam.out.sorted_bam,convert_sam_to_sorted_bam.out.paired,convert_sam_to_sorted_bam.out.sample_name)
    get_mapped_reads(convert_sam_to_sorted_bam.out.sorted_bam,convert_sam_to_sorted_bam.out.paired,convert_sam_to_sorted_bam.out.sample_name)

    emit:
        qc_reads=get_unmapped_reads.out.unmapped_reads
        sample_name=get_unmapped_reads.out.sample_name
        paired=convert_sam_to_sorted_bam.out.paired
}

workflow assembly {
    take:
    sample_name
    reads

    main:
    assemble_with_megahit(sample_name, reads)

    map_reads_fasta_pairs(   
                             assemble_with_megahit.out.sample_name,
                             assemble_with_megahit.out.reads,
                             assemble_with_megahit.out.contigs,
                             assemble_with_megahit.out.paired
                        )
    
    emit:
    contigs=map_reads_fasta_pairs.out.reference_fasta
    sorted_bams=map_reads_fasta_pairs.out.sorted_bam
    sample_name=map_reads_fasta_pairs.out.sample_name
    reads=map_reads_fasta_pairs.out.reads
}
    
workflow binning {
    take:
    sample_name
    sorted_bams
    assembly
    reads
    main:
    get_coverage_for_metabat2(sample_name, sorted_bams, assembly)
    binning_with_metabat2(get_coverage_for_metabat2.out.sample_name,get_coverage_for_metabat2.out.contig,get_coverage_for_metabat2.out.coverage)

    emit:
    metabat2_bins=binning_with_metabat2.out.metabat2_bins
    sample_name=binning_with_metabat2.out.sample_name
    reads=reads

}






// Include modules

include {read_qc_fastp} from "./modules/qc"

include {
    index_bowtie2;
    align_bowtie2;
    convert_sam_to_bam;
    sort_bam;
    convert_sam_to_sorted_bam;
    get_unmapped_reads;
    get_mapped_reads;
    map_reads_fasta_pairs
        } from "./modules/alignment"

include {assemble_with_megahit} from "./modules/assembly"

include {get_coverage_for_metabat2;
        binning_with_metabat2;
        add_prefix_to_fasta} from "./modules/binning"

include {tableToDict;
        concatenate_files;
        } from "./modules/files"

include {get_sequences_from_sra} from "./modules/download"

include {estimate_abundance_coverm} from './modules/abundance'

include {compare_instrain_profiles;
         profile_with_instrain;
         make_stb_file_instrain} from './modules/strain'

include { dereplicate_drep;write_genome_list } from './modules/dereplication'

include {find_genes_prodigal} from './modules/genes'