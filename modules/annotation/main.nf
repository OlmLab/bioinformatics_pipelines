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
    path "gtdbtk_data", emit: gtdbtk_db
    
    script:
    """
    wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
    tar -xvzf gtdbtk_data.tar.gz
    rm gtdbtk_data.tar.gz
    """
}