params.output_dir = "./output"
params.paired = false
params.host_indexed = false
params.binning_extension = "fa"
params.add_fasta_prefix = false
params.is_genome_db=null // default is null
params.is_stb_db=null // default is null
params.is_genes=null // default is null
params.roadmap_5_pairmode="paired"
params.metaphlan_b_distance="bray-curtis"
params.metaphlan_diversity="beta"
params.sylph_db = null
params.sylph_db_link="http://faust.compbio.cs.cmu.edu/sylph-stuff/gtdb-r220-c200-dbv1.syldb"
params.include_metaphlan=false
params.metaphlan_db = null // default is null
params.is_strain_pop_treshold=99
params.is_strain_cos_treshold=99
params.is_strain_con_treshold=99
params.kraken2_db = null // default is null
params.kraken2_kmer_size = 100
params.kraken2_classification_level="S"
params.kraken2_abundance_threshold=10
params.include_kraken2=false
params.kraken2_db_link="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20250402.tar.gz"
params.humann_uniref90 = null // default is null
params.humann_chocophlan = null // default is null
// EXCLUDE PARAMETERS
params.exclude_kraken=false
params.exclude_metaphlan=false
params.exclude_sylph=false
params.exclude_humann=false
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
            if (params.input_reads)
            {
                table=tableToDict(file("${params.input_reads}"))
                reads_1=Channel.fromPath(table["reads1"].collect{t->file(t)})
                reads_2=Channel.fromPath(table["reads2"].collect{t->file(t)})
                inputs=reads_1.merge(reads_2)
                sample_names=Channel.from(table["sample_name"])
                params.roadmap_2_input_type="reads"
            }
            else if (params.input_bams)
            {
                table=tableToDict(file("${params.input_bams}"))
                inputs=Channel.fromPath(table["bam_files"].collect{t->file(t)})
                sample_names=Channel.from(table["sample_name"])
                params.roadmap_2_input_type="bams"
            }
            else
            {
                error "Please provide the reads or bam files information using the input_reads or input_bam parameter."
            }
            
            
            
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
            roadmap_2(sample_names, inputs, fasta_file)
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
            params.roadmap_2_input_type="reads"
            if (params.input_type=="sra")
            {
                table=tableToDict(file("${params.input_file}"))
                get_sequences_from_sra(Channel.fromList(table["Run"]))
                sample_names=get_sequences_from_sra.out.sra_ids
                reads=get_sequences_from_sra.out.fastq_files
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
            if (params.force_genomes)
            {
              force_genomes=Channel.fromPath(file(params.force_genomes))  
            }
            else
            {
              force_genomes=Channel.empty()
            }
            params.roadmap2_input_type="reads"

            roadmap_3_2(sample_names, reads, fasta_file, force_genomes)
    
    }
    else if (params.roadmap_id=="roadmap_5")
    {
        if (!params.input_fastas)
        {
            error "Please provide a genome for mapping reads."
        }
        if (!params.input_reads)
        {
            error "Please provide the reads information using the input_reads parameter."
        }
        table=tableToDict(file("${params.input_reads}"))
        reads_1=Channel.fromPath(table["reads1"].collect{t->file(t)})
        reads_2=Channel.fromPath(table["reads2"].collect{t->file(t)})
        reads=reads_1.merge(reads_2)
        sample_name=Channel.fromList(table["sample_name"])
        samples_reads=sample_name.merge(reads)
        genomes=Channel.fromPath(tableToDict(file("${params.input_fastas}"))["fasta_files"].collect{t->file(t)})
        if (params.roadmap_5_pairmode=="cross")
        {
            samples_reads.combine(genomes).set{inputs}
        }
        else if (params.roadmap_5_pairmode=="paired")
        {
            samples_reads.merge(genomes).set{inputs}
        }
        inputs.multiMap{v->
            sn:v[0]
            rd:[v[1],v[2]]
            gn:v[3]
            tr:true
        }.set{ins}
        roadmap_5(ins.sn, ins.rd, ins.gn,ins.tr)

    }
    else if (params.roadmap_id=="roadmap_6")
    {
        if (params.input_type=="sra")
        {
            table=tableToDict(file("${params.input_file}"))
            get_sequences_from_sra(Channel.fromList(table["Run"]))
            sample_names=get_sequences_from_sra.out.sra_ids
            reads=get_sequences_from_sra.out.fastq_files
            roadmap_6(sample_names, reads)
        }
        else if (params.input_type=="local")
        {
            table=tableToDict(file("${params.input_file}"))
            reads_1=Channel.fromPath(table["reads1"].collect{t->file(t)})
            reads_2=Channel.fromPath(table["reads2"].collect{t->file(t)})
            reads=reads_1.merge(reads_2)
            sample_name=Channel.fromList(table["sample_name"])
            roadmap_6(sample_name, reads)
        }
        else
        {
            error "Please provide the reads information using the file parameter."
        }
    }
    else if (params.roadmap_id=="roadmap_7")
    {
        if (params.bins_dir)
        {
            bins=Channel.fromPath(params.bins_dir,type: 'file')
        }
        else if(params.input_bins_table)
        {
            table=tableToDict(file("${params.input_bins_table}"))
            bins=Channel.fromPath(table["fasta_files"].collect{t->file(t)})

        }
        else
        {
            error "Please provide the bins information using either the bins_dir or input_bins_table parameter."
        }
        roadmap_7(bins.collect())
        
    }
    else if (params.roadmap_id=="roadmap_8")
    {
        if (params.input_type=="sra")
        {
            table=tableToDict(file("${params.input_file}"))
            get_sequences_from_sra(Channel.fromList(table["Run"]))
            sample_names=get_sequences_from_sra.out.sra_ids
            reads=get_sequences_from_sra.out.fastq_files
            host_genome=file(params.host_genome)
            host_genome_gtf=file(params.host_genome_gtf)
            bulk_rna_seq(sample_names, reads, host_genome, host_genome_gtf)
        }
        else if (params.input_type=="local")
        {
            table=tableToDict(file("${params.input_file}"))
            reads_1=Channel.fromPath(table["reads1"].collect{t->file(t)})
            reads_2=Channel.fromPath(table["reads2"].collect{t->file(t)})
            reads=reads_1.merge(reads_2)
            sample_name=Channel.fromList(table["sample_name"])
            host_genome=file(params.host_genome)
            host_genome_gtf=file(params.host_genome_gtf)
            bulk_rna_seq(sample_name, reads, host_genome, host_genome_gtf)
        }
        else
        {
            error "Please provide the reads information using the file parameter."
        }
    }
    else if (params.roadmap_id=="roadmap_9")
    {
        if (params.input_type=="sra")
        {
            table=tableToDict(file("${params.input_file}"))
            get_sequences_from_sra(Channel.fromList(table["Run"]))
            sample_names=get_sequences_from_sra.out.sra_ids
            reads=get_sequences_from_sra.out.fastq_files
        }
        else if (params.input_type=="local")
        {
            table=tableToDict(file("${params.input_file}"))
            reads_1=Channel.fromPath(table["reads1"].collect{t->file(t)})
            reads_2=Channel.fromPath(table["reads2"].collect{t->file(t)})
            reads=reads_1.merge(reads_2)
            sample_name=Channel.fromList(table["sample_name"])
            
        }
        else
        {
            error "Please provide the reads information using the file parameter."
        }
        assemble_rna_spades(sample_name, reads)
        get_circular_contigs_cirit(assemble_rna_spades.out.sample_name, assemble_rna_spades.out.soft_filtered_transcripts)
        
    }

    else if (params.roadmap_id=="download_samples")
    {

        table=tableToDict(file("${params.input_file}"))
        get_sequences_from_sra(Channel.fromList(table["Run"]))

    }
    
    else if (params.roadmap_id=="roadmap_dev")
    {


        table=tableToDict(params.input_profiles)
        profile1=Channel.fromPath(table["profile1"].collect{t->file(t)})
        profile2=Channel.fromPath(table["profile2"].collect{t->file(t)})
        profile_pairs=profile1.merge(profile2)
        stb_file=file(params.stb_file)
        test_customized_compared(profile_pairs,stb_file)

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
    inputs
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
        find_genes_prodigal(genome_db)
        genes=find_genes_prodigal.out.genes_fna
    }
    else
    {
        genes=file(params.is_genes)
    }

    if (params.roadmap_2_input_type=="reads")
    {
        index_bowtie2(genome_db, "genomes_db")
        bowtie2_to_sorted_bam(sample_names, index_bowtie2.out.reference_genome, inputs, index_bowtie2.out.bowtie2_index_files,false)
        sorted_bam=bowtie2_to_sorted_bam.out.sorted_bam
    }
    else if (params.roadmap_2_input_type=="bams")
    {
        sorted_bam=inputs
    }
    profile_with_instrain(sorted_bam,genome_db, stb, genes)
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
    external_genomes
    main:
    roadmap_3(fasta_file.collect())
    dereplicated_genomes=roadmap_3.out.dereplicated_genomes
    dereplicate_drep=dereplicated_genomes.mix(external_genomes).flatten().unique().collect()
    roadmap_2(sample_names, reads, dereplicate_drep)

    emit:
    instrain_profiles=roadmap_2.out.instrain_profiles

}

workflow roadmap_5 {
    //This roadmap is a minimal roadmap for mapping reads a set of reads to a set of genomes 
    take:
    sample_name
    reads
    genome
    paired
    
    main:
    map_reads_fasta_pairs(sample_name, reads, genome,paired)
    get_mapped_reads(map_reads_fasta_pairs.out.sorted_bam,map_reads_fasta_pairs.out.paired,map_reads_fasta_pairs.out.sample_name)
    get_unmapped_reads(map_reads_fasta_pairs.out.sorted_bam,map_reads_fasta_pairs.out.paired,map_reads_fasta_pairs.out.sample_name)

}

workflow roadmap_6{
    // This roadmap analyzes QCed metagenomics reads and performs functional and taxonomic profiling using
    // Reference-based methods.
    take:
    sample_name
    reads
    main:

    sample_name.multiMap{t->
        sample_name_sylph:t
        sample_name_metaphlan:t
    }.set{sample_names}


    reads.multiMap{t->
        reads_sylph:t
        reads_metaphlan:t
    }.set{reads_all}

    reads_all.reads_sylph.multiMap{t->
    reads1_sylph: t[0]
    reads2_sylph: t[1]
    }.set{sylph_reads}

    if (!params.exclude_sylph)
    {
        if (params.sylph_db)
    {
        sylph_db=file(params.sylph_db)
    }
    else
    {
        sylph_db=download_sylph_db()
    }

    estimate_abundance_sylph(sylph_reads.reads1_sylph.collect(), sylph_reads.reads2_sylph.collect(), sylph_db)
    }
   
    if (!params.exclude_metaphlan)
    {
        if (params.metaphlan_db)
          {
           metaphlan_db=file(params.metaphlan_db)
          }
        else
          {
          metaphlan_db=download_metaphlan_db()
          }
        estimate_abundance_metaphlan(sample_names.sample_name_metaphlan,reads_all.reads_metaphlan, metaphlan_db)
        merge_metaphlan_tables(estimate_abundance_metaphlan.out.abundance.collect())
        calculate_diversity_metaphlan(merge_metaphlan_tables.out.merged_table)
    }
   
    if (!params.exclude_kraken)
    {
        if (params.kraken2_db)
        {
            kraken2_db=file(params.kraken2_db)
        }
        else
        {
            kraken2_db=download_kraken2_db(params.kraken2_db_link)
        }
        classify_kraken2(sample_name,reads, kraken2_db)
        estimate_abundance_bracken(classify_kraken2.out.sample_name,classify_kraken2.out.kraken_report,kraken2_db)
    }

    
    if (!params.exclude_humann)
    {
        if (params.humann_uniref90)
        {
            humann_uniref90=file(params.humann_uniref90)
        }
        else
        {
	    download_humann_uniref90()
            humann_uniref90=download_humann_uniref90.out.humann_uniref90
        }
        if (params.humann_chocophlan)
        {
            humann_chocophlan=file(params.humann_chocophlan)
        }
        else
        {
            download_humann_chocophlan()
	    humann_chocophlan=download_humann_chocophlan.out.humann_chocophlan
        }
        profile_humann(sample_name, reads, humann_chocophlan, humann_uniref90, metaphlan_db)    
    }
}

workflow roadmap_7{
    // This roadmap takes a set of bins and performs taxonomic assignment using GTDB-Tk.//
    take:
    bins
    main:
    if (params.gtdbtk_db)
    {
        gtdbtk_db=file(params.gtdbtk_db)
    }
    else
    {
        gtdbtk_db=download_gtdbtk_db()
    }
    assign_taxonomy_gtdb_tk(bins, gtdbtk_db)
    

}

workflow bulk_rna_seq{
    // Standard bulk RNA-Seq workflow for quality control, alignment, and quantification.

    take:
    sample_name
    reads
    host_genome
    host_genome_gtf
    main:
    read_qc_fastp(sample_name, reads)
    index_star(host_genome, host_genome_gtf)
    align_star(sample_name,index_star.out.star_index_files, read_qc_fastp.out.fastp_qcd_reads)

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
    bowtie2_to_sorted_bam(read_qc_fastp.out.sample_name,index_bowtie2.out.reference_genome, read_qc_fastp.out.fastp_qcd_reads,index_bowtie2.out.bowtie2_index_files,true)  
    get_unmapped_reads(bowtie2_to_sorted_bam.out.sorted_bam,bowtie2_to_sorted_bam.out.paired,bowtie2_to_sorted_bam.out.sample_name)
    get_mapped_reads(bowtie2_to_sorted_bam.out.sorted_bam,bowtie2_to_sorted_bam.out.paired,bowtie2_to_sorted_bam.out.sample_name)

    emit:
        qc_reads=get_unmapped_reads.out.unmapped_reads
        sample_name=get_unmapped_reads.out.sample_name
        paired=bowtie2_to_sorted_bam.out.paired
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

workflow test_customized_compared{
    take:
    profile_pairs
    stb_file
    main:
    compare_general_customized(profile_pairs,stb_file)

}




// Include modules

include {read_qc_fastp} from "./modules/qc"

include {
    index_bowtie2;
    align_bowtie2;
    convert_sam_to_bam;
    sort_bam;
    get_unmapped_reads;
    get_mapped_reads;
    map_reads_fasta_pairs;
    bowtie2_to_sorted_bam;
    index_star;
    align_star;
        } from "./modules/alignment"

include {assemble_with_megahit;
        assemble_rna_spades;
        get_circular_contigs_cirit;
        } from "./modules/assembly"

include {get_coverage_for_metabat2;
        binning_with_metabat2;
        add_prefix_to_fasta} from "./modules/binning"

include {tableToDict;
        concatenate_files;
        } from "./modules/files"

include {get_sequences_from_sra} from "./modules/download"

include {estimate_abundance_coverm;
         estimate_abundance_metaphlan;
            estimate_abundance_sylph;
            merge_metaphlan_tables;
            calculate_diversity_metaphlan;
            download_sylph_db;
            download_metaphlan_db;
            download_kraken2_db;
            classify_kraken2;
            estimate_abundance_bracken;
         } from './modules/abundance'

include {compare_instrain_profiles;
         profile_with_instrain;
         make_stb_file_instrain;
         compare_general_customized;
         sample_pairs;
         } from './modules/strain'

include { dereplicate_drep;write_genome_list } from './modules/dereplication'

include {find_genes_prodigal} from './modules/genes'

include {
    download_humann_chocophlan;
    download_humann_uniref90;
    profile_humann;
} from './modules/metabolism'

include {assign_taxonomy_gtdb_tk;
        download_gtdbtk_db;
        } from './modules/annotation'