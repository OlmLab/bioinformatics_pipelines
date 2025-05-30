process download_files {
    /*
    * This process downloads files from a given URL and saves them to a specified directory.
    * @param url: The URL of the file to download.
    * @param output_dir: The directory where the downloaded file will be saved.
    */
    
    publishDir "${params.output_dir}/downloaded_files", mode: 'copy'
    
    input:
    val url
    val output_dir
    
    output:
    path "${output_dir}/*", emit: downloaded_files

    script:
    """
    wget -P ${output_dir} ${url}
    """
}

process get_sequences_from_sra {
    /*
    * This process retrieves sequences from the SRA database using the fastq-dump tool.
    * It takes in a list of SRA IDs and outputs the corresponding FASTQ files.
    * @param sra_ids: SRA ID to retrieve.
    */
    
    publishDir "${params.output_dir}/sra_sequences", mode: 'copy'
    
    input:
    val sra_ids
    
    output:
    path "${sra_ids}/${sra_ids}*.fastq.gz",emit: fastq_files
    val sra_ids, emit: sra_ids
    
    script:
    """
    prefetch ${sra_ids}
    fasterq-dump --split-files --outdir ${sra_ids} ${sra_ids}
    gzip ${sra_ids}/${sra_ids}*.fastq
    rm -rf ${sra_ids}/${sra_ids}.sra
    """
}


