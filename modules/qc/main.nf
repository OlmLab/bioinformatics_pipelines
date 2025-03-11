process read_qc_fastp {
    /*
    * This process applies quality control to the input FASTQ files using fastp. It also generates
    * a report in HTML format and a JSON file with the statistics.
    */
    label "mem_0_cpu_0_time_0"
    publishDir "${params.output_dir}/fastp_qc/${sample_name}", mode: 'copy'
    input:
    val sample_name
    path reads
    output:
    path "${sample_name}_fastp.html", emit: report_html
    path "${sample_name}_fastp.json", emit: report_json
    path("${sample_name}_fastp*.fastq.gz"), emit: fastp_qcd_reads
    val sample_name, emit:sample_name
    script:
    if (reads.size() == 2) {
        """
        fastp \
            -i ${reads[0]} \
            -I ${reads[1]} \
            -o ${sample_name}_fastp_1.fastq.gz \
            -O ${sample_name}_fastp_2.fastq.gz \
            --html ${sample_name}_fastp.html \
            --json ${sample_name}_fastp.json \
            --thread ${task.cpus} \
    """}
    else {
        """
        fastp \
            -i ${reads[0]} \
            -o ${sample_name}_fastp.fastq.gz \
            --html ${sample_name}_fastp.html \
            --json ${sample_name}_fastp.json \
            --thread ${task.cpus} \
        """
    } 

}

  