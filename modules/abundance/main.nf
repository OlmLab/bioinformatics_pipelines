process  estimate_abundance_coverm{
    /*
    * This process estimates the abundance of bins using CoverM.
    * It takes in the bins and reads, and outputs the abundance estimates.
    */
    publishDir "${params.output_dir}/coverm_abundance/${sample_name}", mode: 'copy'
    input:
    val sample_name
    path bins
    path reads
    val extension
    output:
    path "${sample_name}_relative_abundance.tsv", emit: abundance
    script:
    if (reads.size() == 2) {
    """
    coverm genome \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        --genome-fasta-directory . \\
         --genome-fasta-extension ${extension} \\
        --threads ${task.cpus} \\
        --output-file ${sample_name}_relative_abundance.tsv \\
    
    """
    }
    else{
        """
    coverm genome \\
        --single ${reads[0]} \\
        --genome-fasta-directory . \\
         --genome-fasta-extension ${extension} \\
        --threads ${task.cpus} \\
        --output-file ${sample_name}_relative_abundance.tsv \\
\    """
    }
}

process estimate_abundance_sylph{
    /*
    * This process estimates the abundance of bins using Sylph. It takse the reads and prepared database.
    */
    publishDir "${params.output_dir}/sylph_abundance/", mode: 'copy'
    input:
    path reads_1
    path reads_2
    path sylph_db

    output:
    path "sylph_abundance.tsv", emit: abundance
    script:
    """
    sylph profile ${sylph_db} -1 ${reads_1} -2 ${reads_2} -t ${task.cpus} > sylph_abundance.tsv
    """
}

process estimate_abundance_metaphlan{
    /*
    * This process estimates the abundance of bins using MetaPhlAn. It takes in the reads and outputs the abundance estimates.
    */
    publishDir "${params.output_dir}/metaphlan_abundance/", mode: 'copy'
    input:
    val sample_name
    path reads
    path metaphlan_db

    output:
    path "${sample_name}_metaphlan.tsv", emit: abundance
    script:
    """
    metaphlan ${reads[0]},${reads[1]}  --nproc ${task.cpus} --bowtie2db ${metaphlan_db} --input_type fastq -o ${sample_name}_metaphlan.tsv
    """
}

process merge_metaphlan_tables {
    input:
    path abundance
    output:
    path "merged_metaphlan.tsv" ,emit:merged_table
    script:
    """
    merge_metaphlan_tables.py ${abundance} > merged_metaphlan.tsv
    """
}

process calculate_diversity_metaphlan{
    input:
    path table
    output:
    path "*", emit:diversity
    script:
    """
    calculate_diversity.R -f ${table} -d ${params.metaphlan_diversity} -m ${params.metaphlan_b_distance} -o .
    """
}

process download_sylph_db{
    input:
    path sylph_db
    output:
    path "*.syldb", emit: sylph_db
    script:
    """
    wget ${params.sylph_db_link}
    
    """
}
process download_metaphlan_db {
    output:
    path "*", emit: metaphlan_db
    script:
    """
    metaphlan --install --bowtie2db .
    
    """
}