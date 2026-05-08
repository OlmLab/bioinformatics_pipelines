process profile_with_instrain{
    /*
    * This process profiles the bins using InStrain.
    * It takes in the bins and reads, and outputs the profiling results.
    */
    publishDir "${params.output_dir}/instrain/profile/${bamfile.baseName}", mode: params.publish_dir_mode
    input:
    path bamfile
    path fastafile 
    path stb_file
    path genes

    output:
    path("${bamfile.baseName}_instrain_profile"),emit: instrain_profiles
    script:
    def extra_args = task.ext.args ?: ''
    """
    inStrain profile ${bamfile} ${fastafile} -o ${bamfile.baseName}_instrain_profile -p ${task.cpus} -s ${stb_file} --database_mode --gene_file ${genes} ${extra_args}
    """
}

process compare_instrain_profiles{
    /*
    * This process compares the InStrain profiles of the bins.
    * It takes in the profiles and outputs the comparison results.
    */
    publishDir "${params.output_dir}/instrain/compare", mode: params.publish_dir_mode
    input:
    path instrain_profiles 
    path stb_file
    output:
    path "is_compare", emit: instrain_compare
    
    script:
    """
    inStrain compare -i ${instrain_profiles} -o is_compare -p  ${task.cpus} -s ${stb_file} --database_mode

    """
}

process make_stb_file_instrain{
    /*
    * This process creates a stb file for InStrain.
    * It takes in the comparison results and outputs the STB file.
    */
    publishDir "${params.output_dir}/instrain/stb", mode: params.publish_dir_mode
    input:
    path fastafiles
    val name
    output:
    path "${name}.stb", emit: stb_file
    script:
    """
    parse_stb.py --reverse -f ${fastafiles} -o ${name}.stb
    """
}


