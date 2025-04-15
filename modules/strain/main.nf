process profile_with_instrain{
    /*
    * This process profiles the bins using InStrain.
    * It takes in the bins and reads, and outputs the profiling results.
    */
    publishDir "${params.output_dir}/instrain/profile/${bamfile.baseName}", mode: 'copy'
    input:
    path bamfile
    path fastafile 
    path stb_file
    path genes
    output:
    path("${bamfile.baseName}_instrain_profile"),emit: instrain_profiles
    script:
    """
    inStrain profile ${bamfile} ${fastafile} -o ${bamfile.baseName}_instrain_profile -p ${task.cpus} -s ${stb_file} --database_mode --gene_file ${genes}
    """
}

process compare_instrain_profiles{
    /*
    * This process compares the InStrain profiles of the bins.
    * It takes in the profiles and outputs the comparison results.
    */
    publishDir "${params.output_dir}/instrain/compare", mode: 'copy'
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
    publishDir "${params.output_dir}/instrain/stb", mode: 'copy'
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

process compare_general_customized{
    /*
    * This process compares the InStrain profiles using a customized method.
    * It takes in the profiles and outputs the comparison results.
    */
    publishDir "${params.output_dir}/instrain/compare_customized", mode: 'copy'
    input:
    path instrain_profiles 
    path stb_file
    output:
    path "${instrain_profiles[0].baseName}_${instrain_profiles[1].baseName}_compare.json", emit: compare
    
    script:
    """
    compare_cosani.py compare --profile_1 ${instrain_profiles[0]} --profile_2 ${instrain_profiles[1]} --output_file ${instrain_profiles[0].baseName}_${instrain_profiles[1].baseName}_compare.json --stb_file ${stb_file}

    """
}

process get_customized_compared_comps{
    /*
    * This process gets the customized comparison components.
    * It takes in the comparison results and outputs the components.
    */
    publishDir "${params.output_dir}/instrain/compare_customized", mode: 'copy'
    input:
    path comparison
    output:
    path "${comparison.baseName}_strain_share.json", emit: strainshare
    
    script:
    """
    compare_cosani.py stats --compare_profile ${comparison} --output_file ${comparison.baseName}_strain_share.json --strain_pop_treshold ${params.is_strain_pop_treshold} --strain_cos_treshold ${params.is_strain_cos_treshold} --strain_con_treshold ${params.is_strain_con_treshold} 

    """
}