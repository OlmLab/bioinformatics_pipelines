process simulate_short_reads_art_illumina {
    /*
    * Simulate Illumina short reads from a genome using ART (art_illumina).
    * ART uses empirical, platform-specific error profiles (--art_illumina_platform).
    * @param sample_name : sample identifier used as output prefix
    * @param genome      : genome FASTA file (plain or gzipped)
    * Relevant parameters:
    *   --synthetic_paired, --synthetic_read_length, --synthetic_coverage,
    *   --synthetic_fragment_mean, --synthetic_fragment_sd, --synthetic_seed,
    *   --art_illumina_platform, --art_illumina_args
    */
    publishDir "${params.output_dir}/synthetic_reads/${sample_name}", mode: params.publish_dir_mode, pattern: "*.fastq.gz"
    publishDir "${params.output_dir}/synthetic_reads/csv",            mode: params.publish_dir_mode, pattern: "*.csv"

    input:
    val  sample_name
    path genome

    output:
    path "${sample_name}*.fastq.gz", emit: synthetic_reads
    path "${sample_name}.csv",       emit: sample_csv
    val  sample_name,                emit: sample_name

    script:
    def extra_args = task.ext.args ?: ''
    def fastq_dir = "${params.output_dir}/synthetic_reads/${sample_name}"
    def decompress = genome.name.endsWith('.gz') ? "gzip -cd ${genome} > reference.fa" : "ln -s ${genome} reference.fa"
    if (params.synthetic_paired) {
        """
        ${decompress}
        art_illumina \
            -ss ${params.art_illumina_platform} \
            -i reference.fa \
            -p \
            -l ${params.synthetic_read_length} \
            -f ${params.synthetic_coverage} \
            -m ${params.synthetic_fragment_mean} \
            -s ${params.synthetic_fragment_sd} \
            -rs ${params.synthetic_seed} \
            -na \
            ${extra_args} \
            -o ${sample_name}_
        mv ${sample_name}_1.fq ${sample_name}_1.fastq
        mv ${sample_name}_2.fq ${sample_name}_2.fastq
        pigz -p ${task.cpus} ${sample_name}_1.fastq ${sample_name}_2.fastq
        printf 'sample_name,reads1,reads2\\n' > ${sample_name}.csv
        printf '${sample_name},${fastq_dir}/${sample_name}_1.fastq.gz,${fastq_dir}/${sample_name}_2.fastq.gz\\n' >> ${sample_name}.csv
        """
    } else {
        """
        ${decompress}
        art_illumina \
            -ss ${params.art_illumina_platform} \
            -i reference.fa \
            -l ${params.synthetic_read_length} \
            -f ${params.synthetic_coverage} \
            -rs ${params.synthetic_seed} \
            -na \
            ${extra_args} \
            -o ${sample_name}
        mv ${sample_name}.fq ${sample_name}.fastq
        pigz -p ${task.cpus} ${sample_name}.fastq
        printf 'sample_name,reads1\\n' > ${sample_name}.csv
        printf '${sample_name},${fastq_dir}/${sample_name}.fastq.gz\\n' >> ${sample_name}.csv
        """
    }
}

process simulate_short_reads_wgsim {
    /*
    * Simulate short reads from a genome using wgsim.
    * wgsim uses a uniform base error rate and can plant SNPs/indels; the planted
    * mutations are written to <sample_name>_wgsim_mutations.txt.
    * wgsim takes a number of read pairs rather than a coverage, so the number of
    * reads is derived from the genome length and --synthetic_coverage.
    * @param sample_name : sample identifier used as output prefix
    * @param genome      : genome FASTA file (plain or gzipped)
    * Relevant parameters:
    *   --synthetic_paired, --synthetic_read_length, --synthetic_coverage,
    *   --synthetic_fragment_mean, --synthetic_fragment_sd, --synthetic_seed,
    *   --wgsim_error_rate, --wgsim_mutation_rate, --wgsim_indel_fraction, --wgsim_args
    */
    publishDir "${params.output_dir}/synthetic_reads/${sample_name}", mode: params.publish_dir_mode, pattern: "*.{fastq.gz,txt}"
    publishDir "${params.output_dir}/synthetic_reads/csv",            mode: params.publish_dir_mode, pattern: "*.csv"

    input:
    val  sample_name
    path genome

    output:
    path "${sample_name}*.fastq.gz",            emit: synthetic_reads
    path "${sample_name}_wgsim_mutations.txt",  emit: mutations
    path "${sample_name}.csv",                  emit: sample_csv
    val  sample_name,                           emit: sample_name

    script:
    def extra_args = task.ext.args ?: ''
    def fastq_dir = "${params.output_dir}/synthetic_reads/${sample_name}"
    def decompress = genome.name.endsWith('.gz') ? "gzip -cd ${genome} > reference.fa" : "ln -s ${genome} reference.fa"
    def reads_per_fragment = params.synthetic_paired ? 2 : 1
    def wgsim_common = """-e ${params.wgsim_error_rate} \
            -r ${params.wgsim_mutation_rate} \
            -R ${params.wgsim_indel_fraction} \
            -1 ${params.synthetic_read_length} \
            -2 ${params.synthetic_read_length} \
            -d ${params.synthetic_fragment_mean} \
            -s ${params.synthetic_fragment_sd} \
            -S ${params.synthetic_seed} \
            ${extra_args}"""
    def count_reads = """
        genome_length=\$(grep -v '^>' reference.fa | tr -d '\\n\\r' | wc -c)
        n_reads=\$(( (genome_length * ${params.synthetic_coverage} + ${reads_per_fragment * params.synthetic_read_length} - 1) / ${reads_per_fragment * params.synthetic_read_length} ))
        """
    if (params.synthetic_paired) {
        """
        ${decompress}
        ${count_reads}
        wgsim \
            -N \$n_reads \
            ${wgsim_common} \
            reference.fa ${sample_name}_1.fastq ${sample_name}_2.fastq > ${sample_name}_wgsim_mutations.txt
        pigz -p ${task.cpus} ${sample_name}_1.fastq ${sample_name}_2.fastq
        printf 'sample_name,reads1,reads2\\n' > ${sample_name}.csv
        printf '${sample_name},${fastq_dir}/${sample_name}_1.fastq.gz,${fastq_dir}/${sample_name}_2.fastq.gz\\n' >> ${sample_name}.csv
        """
    } else {
        """
        ${decompress}
        ${count_reads}
        wgsim \
            -N \$n_reads \
            ${wgsim_common} \
            reference.fa ${sample_name}.fastq /dev/null > ${sample_name}_wgsim_mutations.txt
        pigz -p ${task.cpus} ${sample_name}.fastq
        printf 'sample_name,reads1\\n' > ${sample_name}.csv
        printf '${sample_name},${fastq_dir}/${sample_name}.fastq.gz\\n' >> ${sample_name}.csv
        """
    }
}
