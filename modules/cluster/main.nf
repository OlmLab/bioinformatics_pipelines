
process create_mmseqs_db{
    /*
    * This process creates an MMseqs2 database from the input FASTA file.
    * @param input_fasta: The input FASTA file containing sequences.
    */
    publishDir "${params.output_dir}/mmseqs2_db", mode: 'copy'
    
    input:
    path input_fasta

    output:
    path "${input_fasta}.*", emit: seq_db

    script:
    """
    mmseqs createdb ${input_fasta} ${input_fasta}
    """
}

process mmseqs_linclust{
    /*
    * This process clusters sequences using MMseqs2 linclust algorithm.
    * @param input_fasta: The input FASTA file containing sequences to be clustered.
    * @param identity: The sequence identity threshold for clustering.
    */
    publishDir "${params.output_dir}/mmseqs2", mode: 'copy'
    
    input:
    path input_fasta
    path seq_db

    output:
    path "${input_fasta.fileName}_representative_${params.mmseqs_linclust_identity}.fasta", emit: clustered_sequences

    script:
    """
    mmseqs linclust ${seq_db} clustered_db tmp --min-seq-id ${params.mmseqs_linclust_identity} -c ${params.mmseqs_linclust_coverage} --cov-mode 1 --threads ${task.cpus}
    mmseqs createsubdb clustered_db ${seq_db} clustered_db_rep
    mmseqs convert2fasta clustered_db_rep ${input_fasta.fileName}_representative_${params.mmseqs_linclust_identity}.fasta
    """
}