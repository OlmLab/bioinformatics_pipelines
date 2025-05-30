params.drep_s_ani= "0.95"

process  write_genome_list {

    publishDir "${params.output_dir}/dereplication_drep", mode: 'copy'
    input:
    path genomes
    output:
    path "genome_list.txt", emit: genomes_list
    script:
    """
    ls * > genome_list.txt
    """
}
process dereplicate_drep{

    publishDir "${params.output_dir}/dereplication_drep", mode: 'copy'
    
    input:
    path genomes
    path genomes_list

    output:
    path "drep_output/dereplicated_genomes/*", emit: dereplicated_genomes
    path "drep_output/"
    script:
    def extra_weight_arg= params.drep_extra_weight_table  ? "--extra_weight_table ${params.drep_extra_weight_table}" : ""
    """

    dRep dereplicate drep_output/ -p ${task.cpus} -g ${genomes_list}  --S_ani ${params.drep_s_ani} ${extra_weight_arg}

    """
}