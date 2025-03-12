params.output_dir = "./output_is"
params.paired = false
params.host_indexed = false
params.binning_extension = "fa"
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
            reads_1=table["reads1"]
            reads_2=table["reads2"]
            reads=[reads_1,reads_2].transpose()
            sample_names=table["sample_name"].toList()
            fasta_file=tableToDict(file("${params.input_fastas}"))["fasta_files"].toList()
            roadmap_2(sample_names, reads, fasta_file)
        }
        


        
    }
    

    



workflow seqs_from_sra {
    take:
    sra_ids

    main:
    get_sequences_from_sra(sra_ids)

    emit:
    fastq_files=get_sequences_from_sra.out.fastq_files
    sample_name=get_sequences_from_sra.out.sra_ids
}


workflow roadmap_1{
    // This roadmap takes the raw sequencing reads and performs:
    // quality control, assembly, and binning.
    take:
    sample_name
    reads
    host_genome

    main:
    quality_control(sample_name, reads, host_genome)
    assembly(quality_control.out.sample_name, quality_control.out.qc_reads)
    binning(assembly.out.sample_name, assembly.out.sorted_bams, assembly.out.contigs,assembly.out.reads)
    estimate_abundance_coverm(binning.out.sample_name, binning.out.metabat2_bins, binning.out.reads, params.binning_extension)    
    emit:
    metabat2_bins=binning.out.metabat2_bins
}

//
workflow roadmap_2 {
    take:
    sample_names
    reads
    fasta_file
    main:
    fasta_file_ch=Channel.fromList(fasta_file).map{t->file(t)}
    sample_names_ch=Channel.fromList(sample_names)
    reads_ch=Channel.fromList(reads).map{t->tuple(file(t[0]),file(t[1]))}
    index_bowtie2(fasta_file_ch,sample_names_ch)
    samples_reads=sample_names_ch.merge(reads_ch)
    samples_reads.combine(index_bowtie2.out.bowtie2_index_files.merge(index_bowtie2.out.reference_genome)).multiMap{t->
    sample_names_:t[0]
    reads_:tuple(t[1],t[2])
    index_bowtie2_:tuple(t[3],t[4],t[5],t[6],t[7],t[8])
    reference_genome_:tuple(t[-1])
    }.set{split_results}
    align_bowtie2(split_results.sample_names_,split_results.reference_genome_,split_results.reads_,split_results.index_bowtie2_)
    convert_sam_to_sorted_bam(align_bowtie2.out.bowtie2_sam,align_bowtie2.out.sample_name,align_bowtie2.out.paired)
    profile_with_instrain(convert_sam_to_sorted_bam.out.sorted_bam.merge(split_results.reference_genome_))
    profile_with_instrain.out.instrain_profiles.groupTuple(by:0).set{profiles}
    compare_instrain_profiles(profiles)
}



    // compare_instrain_profiles(profiles, "${fasta_file.baseName}_instrain_compare")

    

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
    index_bowtie2(host_genome,sample_name)
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

include {get_coverage_for_metabat2;binning_with_metabat2} from "./modules/binning"

include {tableToDict} from "./modules/files"

include {get_sequences_from_sra} from "./modules/download"

include {estimate_abundance_coverm} from './modules/abundance'

include {compare_instrain_profiles;profile_with_instrain} from './modules/strain'
