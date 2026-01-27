process assign_taxonomy_gtdb_tk {
    /*
    * This process assigns taxonomy to the input bins using GTDB-Tk.
    */
    publishDir "${params.output_dir}/gtdbtk", mode: 'copy'
    input:
    path bins
    env "GTDBTK_DATA_PATH" // GTDB-Tk database path

    output:
    path "gtdb_annotation", emit: gtdb_annotation
    path "mash_db", emit: mash_db
    script:
    """
    mv ${bins} bins
    gtdbtk classify_wf --genome_dir bins --out_dir gtdb_annotation  --cpus ${task.cpus} -x ${params.binning_extension} --mash_db mash_db

    """
}

process download_gtdbtk_db {
    /*
    * This process downloads the GTDB-Tk database from a given URL and extracts it.
    * @param url: The URL of the GTDB-Tk database to download.
    */
    
    publishDir "${params.output_dir}/gtdbtk_db"
    output:
    path "release*", emit: gtdbtk_db
    
    script:
    """
    wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
    tar -xvzf gtdbtk_data.tar.gz
    rm gtdbtk_data.tar.gz
    """
}

process classify_kraken2_contigs{ 
    publishDir "${params.output_dir}/kraken2/${sample_name}", mode: 'copy'
    input:
    val sample_name
    path contigs_fasta
    path kraken2_db
    output:
    val sample_name, emit: sample_name
    path "${sample_name}_kraken2_*.txt", emit: kraken_report
    script:
    {
    """
    kraken2 --db ${kraken2_db} --threads ${task.cpus} -input ${contigs_fasta} --report ${sample_name}_kraken2_report.txt --output ${sample_name}_kraken2_output.txt  
    """
    }
 
}

process download_eggnog_db {
    publishDir "${params.output_dir}/eggnog_db", mode: 'link'
    output:
    path "eggnog_data", emit: eggnog_db
    script:
    """
    mkdir eggnog_data
    download_eggnog_data.py --data_dir eggnog_data -H -d ${params.eggnog_db_taxonomic_scope} -y
    """
}

process eggnog_annotation {
    /*
    * This process annotates the input contigs using eggNOG-mapper.
    */
    publishDir "${params.output_dir}/eggnog_annotation", mode: 'copy'
    input:
    path genes_fasta
    path eggnog_data_dir

    output:
    path "${genes_fasta.baseName}_eggnog_annotation", emit: eggnog_annotation

    script:
    """
    mkdir ${genes_fasta.baseName}_eggnog_annotation
    emapper.py -i ${genes_fasta} --itype CDS --translate  --output_dir ${genes_fasta.baseName}_eggnog_annotation -o ${genes_fasta.baseName}_eggnog_annotation --cpu ${task.cpus} --data_dir ${eggnog_data_dir} -d ${params.eggnog_db_taxonomic_scope}
    """
}


process download_genomad_db {
    /*
    * This process downloads the Genomad database from a given URL and extracts it.
    * @param url: The URL of the Genomad database to download.
    */
    
    publishDir "${params.output_dir}/"
    output:
    path "genomad_db", emit: genomad_db
    
    script:
    """
    genomad download-database .
    """
}

process annotate_contig_genomad{
    /*
    * This process annotates contigs using Genomad.
    */
    publishDir "${params.output_dir}/genomad_annotation"
    input:
    val sample_name
    path contigs_fasta
    path genomad_db
    

    output:
    path "${sample_name}_genomad_annotation", emit: genomad_annotation

    script:
    """
    mkdir ${sample_name}_genomad_annotation
    genomad end-to-end --cleanup --splits 8 ${contigs_fasta} ${sample_name}_genomad_annotation ${genomad_db}"""
}