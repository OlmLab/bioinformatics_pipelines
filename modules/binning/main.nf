process  get_coverage_for_metabat2
{
    /*
    * This process calculates the coverage of the input assembly using the provided reads. It generates
    * a coverage file that will be used for binning in the next step.
    */
    publishDir "${params.output_dir}/binning/metabat2/${sample_name}", mode: 'copy'
    input:
    val sample_name
    path sorted_bam
    path contig
    output:
    path "${sample_name}_coverage.txt", emit: coverage
    val sample_name, emit: sample_name
    path contig, emit: contig
    script:
    """
    jgi_summarize_bam_contig_depths \\
        --outputDepth ${sample_name}_coverage.txt \\
        ${sorted_bam} 
    """
}



process binning_with_metabat2
{
    /*
    * This process performs binning using MetaBAT2. It takes the input assembly and the coverage
    * information from the previous step and generates bins in a specified output directory.
    */
    publishDir "${params.output_dir}/binning/metabat2/${sample_name}", mode: 'copy'
    input:
    val sample_name
    path assembly
    path coverage
    output:
    path  "${sample_name}_metabat2_bins*", emit: metabat2_bins
    val sample_name, emit: sample_name
    script:
    """ 
    metabat2 \\
        -i ${assembly} \\
        -a ${coverage} \\
        --numThreads ${task.cpus} \\
        --outFile ${sample_name}_metabat2_bins \\
    """
}

process add_prefix_to_fasta {
    input:
    path fasta_file
    output:
    path "prefixed_fasta/*", emit: prefixed_fasta
    script:
    """
    add_prefix_to_fasta.py \\
    ${fasta_file} \\
    --output-dir prefixed_fasta 
    """
}