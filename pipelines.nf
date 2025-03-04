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
            roadmap_1(get_sequences_from_sra.sra_ids, get_sequences_from_sra.out.fastq_files, file(params.host_genome))
        }
        if (params.input_type=="local")
        {
            table=tableToDict(file("${params.input_file}"))
            reads_1=Channel.fromPath(table["reads1"])
            reads_2=Channel.fromPath(table["reads2"])
            reads=reads_1.combine(reads_2)
            sample_name=Channel.fromList(table["sample_name"])
            roadmap_1(sample_name, reads, file(params.host_genome))
        }

        
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
    assembly(sample_name, quality_control.out.qc_reads)
    binning(sample_name, assembly.out.sorted_bams, assembly.out.contigs)
    estimate_abundance_coverm(sample_name, binning.out.metabat2_bins, quality_control.out.qc_reads, params.binning_extension)    
    emit:
    metabat2_bins=binning.out.metabat2_bins
}

//

workflow compare_samples_instrain {
    take:
    reads
    sample_name
    fasta_file
    main:
    index_bowtie2(fasta_file)
    align_bowtie2(sample_name,index_bowtie2.out.reference_genome, reads,index_bowtie2.out.bowtie2_index_files)
    convert_sam_to_sorted_bam(align_bowtie2.out.bowtie2_sam)
    profile_with_instrain(convert_sam_to_sorted_bam.out.sorted_bam, fasta_file)
    // compare_instrain_profiles(profiles, "${fasta_file.baseName}_instrain_compare")

    
}
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
    index_bowtie2(host_genome)
    if (params.paired)
    {
        

        if (params.paired)
        {
          align_bowtie2(sample_name,index_bowtie2.out.reference_genome, read_qc_fastp.out.fastp_qcd_reads,index_bowtie2.out.bowtie2_index_files)  
          convert_sam_to_sorted_bam(align_bowtie2.out.bowtie2_sam)
          get_unmapped_reads(convert_sam_to_sorted_bam.out.sorted_bam,align_bowtie2.out.paired)
          get_mapped_reads(convert_sam_to_sorted_bam.out.sorted_bam,align_bowtie2.out.paired)
        } }    
      
    else
        {
            align_bowtie2(sample_name,index_bowtie2.out.reference_genome, read_qc_fastp.out.fastp_qcd_reads,index_bowtie2.out.bowtie2_index_files)  
            convert_sam_to_sorted_bam(align_bowtie2.out.bowtie2_sam)
            get_unmapped_reads(convert_sam_to_sorted_bam.out.sorted_bam,align_bowtie2.out.paired)
            get_mapped_reads(convert_sam_to_sorted_bam.out.sorted_bam,align_bowtie2.out.paired)
        }
        emit:
        qc_reads=get_unmapped_reads.out.unmapped_reads


}

workflow assembly {
    take:
    sample_name
    reads

    main:
    assemble_with_megahit(sample_name, reads)
    index_bowtie2(assemble_with_megahit.out.contigs)
    align_bowtie2(sample_name,index_bowtie2.out.reference_genome, reads,index_bowtie2.out.bowtie2_index_files)
    convert_sam_to_sorted_bam(align_bowtie2.out.bowtie2_sam)
    
    emit:
    contigs=assemble_with_megahit.out.contigs
    sorted_bams=convert_sam_to_sorted_bam.out.sorted_bam
}
    
workflow binning {
    take:
    sample_name
    sorted_bams
    assembly
    main:
    get_coverage_for_metabat2(sample_name, sorted_bams)
    binning_with_metabat2(sample_name,assembly,get_coverage_for_metabat2.out.coverage)

    emit:
    metabat2_bins=binning_with_metabat2.out.metabat2_bins

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
        } from "./modules/alignment"

include {assemble_with_megahit} from "./modules/assembly"

include {get_coverage_for_metabat2;binning_with_metabat2} from "./modules/binning"

include {tableToDict} from "./modules/files"

include {get_sequences_from_sra} from "./modules/download"

include {estimate_abundance_coverm} from './modules/abundance'

include {compare_instrain_profiles;profile_with_instrain} from './modules/strain'
