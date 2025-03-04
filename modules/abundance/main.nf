process  estimate_abundance_coverm{
    /*
    * This process estimates the abundance of bins using CoverM.
    * It takes in the bins and reads, and outputs the abundance estimates.
    */
    publishDir "${params.output_dir}/coverm_abundance/${sample_name}", mode: 'copy'
    input:
    val sample_name
    path bins
    path reads
    val extension
    output:
    path "${sample_name}_abundance.tsv", emit: abundance
    script:
    if (reads.size() == 2) {
    """
    coverm genome \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        --genome-fasta-directory . \\
         --genome-fasta-extension ${extension} \\
        --threads ${task.cpus} \\
        --output-file ${sample_name}_relative_abundance.tsv \\
    
    """
    }
    else{
        """
    coverm genome \\
        --single ${reads[0]} \\
        --genome-fasta-directory . \\
         --genome-fasta-extension ${extension} \\
        --threads ${task.cpus} \\
        --output-file ${sample_name}_relative_abundance.tsv \\
\    """
    }
}