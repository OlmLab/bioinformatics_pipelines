process assemble_with_megahit{
    /*
    * This process assembles reads using MEGAHIT.
    * It takes in paired-end or single-end reads and outputs the assembled contigs.
    * @param read1: The first read file (or single-end read).
    * @param read2: The second read file (optional, for paired-end reads).
    */

    publishDir "${params.output_dir}/megahit", mode: 'copy'
    input:
    val sample_name 
    path reads
    output:
    path "megahit_out/${sample_name}.contigs.fa", emit: contigs
    path "megahit_out/*"
    
    script:
    if (params.paired) {
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