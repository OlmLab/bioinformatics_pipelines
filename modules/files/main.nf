process concatenate_files {
    /*
    * This process takes in a list of files and concatenates them into a single file.
    */
    publishDir params.output_dir, mode: 'copy'
    input:
    path file_list
    val output_file

    output:
    path "${output_file}", emit: concatenated_file

    script:
    """
    cat ${file_list} > ${output_file}
    """
}
process subsample_reads_reformat {
    /*
    * Randomly subsample reads at a given fraction using BBTools reformat.sh.
    * Paired-end reads are processed together so read pairing is guaranteed to
    * be preserved.
    * @param sample_name : sample identifier used as output prefix
    * @param reads       : one (single-end) or two (paired-end) FASTQ files
    * @param fraction    : fraction of reads to retain, between 0 and 1
    *                      (e.g. 0.1 keeps 10 % of reads)
    */
    publishDir "${params.output_dir}/subsampled_reads/${sample_name}/fraction_${fraction}", mode: 'move', pattern: "*.fastq.gz"
    publishDir "${params.output_dir}/subsampled_reads/csv",                                 mode: 'move', pattern: "*.csv"

    input:
    val  sample_name
    path reads
    val  fraction

    output:
    path "${sample_name}_${fraction}*.fastq.gz", emit: subsampled_reads
    path "${sample_name}_${fraction}.csv",       emit: sample_csv
    val  sample_name,                            emit: sample_name
    val  fraction,                               emit: fraction

    script:
    def fastq_dir = "${params.output_dir}/subsampled_reads/${sample_name}/fraction_${fraction}"
    if (reads.size() == 2) {
        """
        reformat.sh \
            in1=${reads[0]} \
            in2=${reads[1]} \
            out1=${sample_name}_${fraction}_1.fastq.gz \
            out2=${sample_name}_${fraction}_2.fastq.gz \
            samplerate=${fraction} \
            sampleseed=${params.subsample_seed}
        printf 'sample_name,reads1,reads2\\n' > ${sample_name}_${fraction}.csv
        printf '${sample_name},${fastq_dir}/${sample_name}_${fraction}_1.fastq.gz,${fastq_dir}/${sample_name}_${fraction}_2.fastq.gz\\n' >> ${sample_name}_${fraction}.csv
        """
    } else {
        """
        reformat.sh \
            in=${reads[0]} \
            out=${sample_name}_${fraction}.fastq.gz \
            samplerate=${fraction} \
            sampleseed=${params.subsample_seed}
        printf 'sample_name,reads1\\n' > ${sample_name}_${fraction}.csv
        printf '${sample_name},${fastq_dir}/${sample_name}_${fraction}.fastq.gz\\n' >> ${sample_name}_${fraction}.csv
        """
    }
}

def tableToDict(file, delimiter = ',') {
    /*
    * This function reads a CSV file and converts it into a dictionary (map) where the keys are the headers
    * and the values are lists of the corresponding column values.
    * @param file: The CSV file to read.
    * @param delimiter: The delimiter used in the CSV file (default is comma).
    * @return: A map where keys are headers and values are lists of column values.
    */



    def result = [:]
    def lines = file.readLines()
    
    // Get all headers
    def headers = lines[0].split(delimiter).collect { it.trim() }
    
    // Initialize empty arrays for each header
    headers.each { header ->
        result[header] = []
    }
    
    // Collect all values for each column
    lines[1..-1].each { line ->
        def values = line.split(delimiter).collect { it.trim() }
        headers.eachWithIndex { header, i ->
            result[header] << values[i]
        }
    }
    
    return result
}
