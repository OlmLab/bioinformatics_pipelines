def write_genome_list(genomes) {
    def genome_list_file = file("${params.output_dir}/dereplication_drep/genome_list.txt")
    (genome_list_file.parent).mkdirs() 
    genome_list_file.text = genomes.collect { genome ->
        genome.getName()
    }.join("\n")
    return genome_list_file
}

process dereplicate_drep{
    publishDir "${params.output_dir}/dereplication_drep"
    input:
    path genomes
    path genomes_list

    output:
    path "drep_output/dereplicated_genomes/*", emit: dereplicated_genomes

    script:

    """
    dRep dereplicate drep_output/ -p ${task.cpus} -g ${genomes_list}  

    """
}