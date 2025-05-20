process find_genes_prodigal{
    publishDir "${params.output_dir}/genes/prodigal", mode: 'copy'
    input:
    path genome
    output:
    path "${genome.name}_genes.fna", emit: genes_fna
    path "${genome.name}_genes.faa", emit: genes_faa
    script:
    """
    prodigal -i ${genome} -d ${genome.name}_genes.fna -a ${genome.name}_genes.faa -p meta
    
    """
}

