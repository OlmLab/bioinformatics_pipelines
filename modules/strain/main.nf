process profile_with_instrain{
    /*
    * This process profiles the bins using InStrain.
    * It takes in the bins and reads, and outputs the profiling results.
    */
    publishDir "${params.output_dir}/instrain/${fasta.baseName}"
    input:
    tuple(path(bamfile),path(fasta))   
    output:
    tuple val(fasta.baseName),path("${bamfile.baseName}_${fasta.baseName}_instrain_profile"),emit: instrain_profiles
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
    publishDir "${params.output_dir}/instrain/${fastafile}"
    input:
    tuple val(fastafile),path(instrain_profiles)
    output:
    path "${fastafile}_compare", emit: instrain_compare
    
    script:
    """
    inStrain compare -i ${instrain_profiles} -o ${fastafile}_compare -p ${task.cpus}
    """
}