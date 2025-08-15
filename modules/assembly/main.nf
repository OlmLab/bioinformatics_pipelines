process assemble_with_megahit{
    /*
    * This process assembles reads using MEGAHIT.
    * It takes in paired-end or single-end reads and outputs the assembled contigs.
    * @param read1: The first read file (or single-end read).
    * @param read2: The second read file (optional, for paired-end reads).
    */

    publishDir "${params.output_dir}/megahit/${sample_name}", mode: 'copy'
    input:
    val sample_name 
    path reads
    output:
    path "megahit_out/${sample_name}.contigs.fa", emit: contigs
    path "megahit_out/*"
    val sample_name, emit: sample_name
    path reads, emit: reads
    val paired , emit: paired
    
    script:
    paired=(reads.size() == 2)
    if (reads.size() == 2) {
        """
        megahit \
            -1 ${reads[0]} \
            -2 ${reads[1]} \
            -o megahit_out \
            -t ${task.cpus} \
            --out-prefix ${sample_name}
        """
    } else {
        """
        megahit \
            -r ${reads[0]} \
            -o megahit_out \
            -t ${task.cpus} \
            --out-prefix ${sample_name}
        """
    }
}

process assemble_rna_spades{
publishDir "rna_spades/${sample_name}", mode: params.publish_mode
input:
    val sample_name
    path reads
output:
    path "hard_filtered_transcripts.fasta", emit: hard_filtered_transcripts
    path "soft_filtered_transcripts.fasta", emit: soft_filtered_transcripts
    val sample_name, emit: sample_name

script:
if (reads.size() == 2) {
    """
     rnaspades.py -1 ${reads[0]} -2 ${reads[1]} -o . --threads ${task.cpus} 
    """
}
else{
    """
     rnaspades.py -s ${reads[0]} -o . --threads ${task.cpus} 
    """
}
}

process get_circular_contigs_cirit{
publishDir "rna_spades/${sample_name}", mode: params.publish_mode
input:
    val sample_name
    path assembly
output:
    path "${sample_name}_circular_contigs.fasta", emit: circular_contigs
script:
    """
    java -jar /Cirit.jar -i ${assembly} -o ${sample_name}_circular_contigs.fasta 
    """
}