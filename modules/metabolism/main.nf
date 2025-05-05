process download_humann_database {
    /*
    * This process downloads the HUMAnN database from a given URL and extracts it.
    * @param url: The URL of the HUMAnN database to download.
    */
    
    publishDir "${params.output_dir}/humann_databases"
    output:
    path "chocoplan", emit: humann_chocophlan
    path "uniref90", emit: humann_uniref90

    
    script:
    """
    humann_databases --download chocophlan full ./chocophlan
    humann_databases --download uniref uniref90_diamond ./uniref90
    """
}

process profile_humann {
    /*
    * This process profiles the HUMAnN database using the input reads.
    * @param reads: The input reads to profile.
    * @param chocoplan: The path to the chocoplan database.
    * @param uniref90: The path to the uniref90 database.
    */
    
    publishDir "${params.output_dir}/humann_profile", mode: 'copy'
    input:
    val sample_name
    path reads
    output:
    path "*.tsv", emit: humann_profile
    
    script:
    """
    cat ${reads} > ${sample_name}.fastq
    humann --input ${sample_name}.fastq --output . --database chocoplan --threads ${task.cpus} --protein-database ${params.humann_uniref90} --nucleotide-database ${params.humann_chocophlan} --output-format tsv
    rm ${sample_name}.fastq
    """
}