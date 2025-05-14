process download_humann_chocophlan {
    /*
    * This process downloads the HUMAnN database from a given URL and extracts it.
    * @param url: The URL of the HUMAnN database to download.
    */
    
    publishDir "${params.output_dir}/humann_databases"
    output:
    path "chocophlan/chocophlan", emit: humann_chocophlan

    
    script:
    """
    humann_databases --download chocophlan full ./chocophlan --update-config no
    """
}

process download_humann_uniref90 {
    /*
    * This process downloads the HUMAnN database from a given URL and extracts it.
    * @param url: The URL of the HUMAnN database to download.
    */
    
    publishDir "${params.output_dir}/humann_databases"
    output:
    path "uniref90/uniref", emit: humann_uniref90

    
    script:
    """
    humann_databases --download uniref uniref90_diamond ./uniref90 --update-config no
    """
}

process profile_humann {
    /*
    * This process profiles the HUMAnN database using the input reads.
    * @param reads: The input reads to profile.
    * @param sample_name: The name of the sample.
    */
    
    publishDir "${params.output_dir}/humann_profile", mode: 'copy'
    input:
    val sample_name
    path reads
    path humann_chocophlan
    path humann_uniref90 
    path metaphlan_db
    output:
    path "*.tsv", emit: humann_profile
    
    script:
    """
    cat ${reads} > ${sample_name}.fastq.gz
    humann --input ${sample_name}.fastq.gz --output .  --threads ${task.cpus} --protein-database ${humann_uniref90} --nucleotide-database ${humann_chocophlan} --output-format tsv --metaphlan-options "--bowtie2db ${metaphlan_db.name}"
    rm ${sample_name}.fastq.gz
    """
}