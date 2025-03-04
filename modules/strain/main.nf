process profile_with_instrain{
    /*
    * This process profiles the bins using InStrain.
    * It takes in the bins and reads, and outputs the profiling results.
    */
    publishDir "${params.output_dir}/instrain/${bamfile.baseName}_${fasta.baseName}"
    input:
    path bamfile
    path fasta
    
    output:
    path "${bamfile.baseName}_${fasta.baseName}_instrain_profile", emit: instrain_profile
    script:
    """
    inStrain profile ${bamfile} ${fasta} -o ${bamfile.baseName}_${fasta.baseName}_instrain_profile -p ${task.cpus}
    """
}

process compare_instrain_profiles{
    /*
    * This process compares the InStrain profiles of the bins.
    * It takes in the profiles and outputs the comparison results.
    */
    publishDir "${params.output_dir}/instrain/compare"
    input:
    path instrain_profiles
    val compare_name
    output:
    path "${compare_name}", emit: instrain_compare
    
    script:
    """
    inStrain compare ${instrain_profiles} -o ${compare_name} -p ${task.cpus}
    """
}