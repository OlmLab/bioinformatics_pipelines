
nextflow.enable.dsl = 2
// Bowtie2 parameters
params.bowtie2_non_competitive_mapping=false
//
params.single_end = false
params.output_dir = "./output"
params.paired = false
params.host_indexed = false
params.binning_extension = "fa"
params.add_fasta_prefix = false
params.is_genome_db=null // default is null
params.is_stb_db=null // default is null
params.is_genes=null // default is null
params.roadmap_5_pairmode="paired"
params.read_type="short"              // Sequencing technology: short, nanopore, pacbio_clr, pacbio_hifi
params.keep_unmapped_reads=false      // Include unmapped reads in output BAM (default: mapped-only BAM)
params.get_mapped_reads=false         // Additionally output a FASTQ of mapped reads
params.get_unmapped_reads=false       // Additionally output a FASTQ of unmapped reads
params.dump_report=true               // Write pipeline_info/trace.txt and pipeline_info/run_stats.txt
params.metaphlan_b_distance="bray-curtis"
params.metaphlan_diversity="beta"
params.sylph_db = null
params.sylph_db_link="http://faust.compbio.cs.cmu.edu/sylph-stuff/gtdb-r220-c200-dbv1.syldb"
params.include_metaphlan=false
params.metaphlan_db = null // default is null
params.is_strain_pop_treshold=99
params.is_strain_cos_treshold=99
params.is_strain_con_treshold=99
params.kraken2_db = null // default is null
params.kraken2_kmer_size = 100
params.kraken2_classification_level="S"
params.kraken2_abundance_threshold=10
params.include_kraken2=false
params.kraken2_db_link="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20250402.tar.gz"
params.humann_uniref90 = null // default is null
params.humann_chocophlan = null // default is null
// EXCLUDE PARAMETERS
params.exclude_kraken=false
params.exclude_metaphlan=false
params.exclude_sylph=false
params.exclude_humann=false
params.scan_genome_batch_size=10
params.subsample_seed=42
params.fractions=null


// mmseqs2 parameters
params.mmseqs_linclust_identity = 0.95
params.mmseqs_linclust_coverage = 0.8
// ###### MAIN WORKFLOW ###### //
params.build_gene_db_mode="nucleotide" 
params.eggnog_db_taxonomic_scope="2"

def parseLongMetric(value) {
    if (value == null) {
        return null
    }
    def text = value.toString().trim()
    if (!text || text.equalsIgnoreCase('NA')) {
        return null
    }
    try {
        return new BigDecimal(text).longValue()
    } catch (Exception ignored) {
        return null
    }
}

def parseDecimalMetric(value) {
    if (value == null) {
        return null
    }
    def text = value.toString().trim()
    if (!text || text.equalsIgnoreCase('NA')) {
        return null
    }
    text = text.replace('%', '')
    try {
        return new BigDecimal(text)
    } catch (Exception ignored) {
        return null
    }
}

def formatDurationMs(value) {
    Long millis = parseLongMetric(value)
    if (millis == null) {
        return "NA"
    }
    if (millis < 1000L) {
        return "${millis} ms"
    }

    long totalSeconds = (long) (millis / 1000L)
    long days = (long) (totalSeconds / 86400L)
    long hours = (long) ((totalSeconds % 86400L) / 3600L)
    long minutes = (long) ((totalSeconds % 3600L) / 60L)
    long seconds = (long) (totalSeconds % 60L)

    def parts = []
    if (days) {
        parts << "${days}d"
    }
    if (hours) {
        parts << "${hours}h"
    }
    if (minutes) {
        parts << "${minutes}m"
    }
    if (seconds || !parts) {
        parts << "${seconds}s"
    }
    return parts.join(' ')
}

def formatBytes(value) {
    Long bytes = parseLongMetric(value)
    if (bytes == null) {
        return "NA"
    }

    def units = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']
    BigDecimal size = new BigDecimal(bytes)
    int unitIndex = 0
    while (size >= 1024 && unitIndex < units.size() - 1) {
        size = size.divide(new BigDecimal(1024), 2, java.math.RoundingMode.HALF_UP)
        unitIndex += 1
    }
    return "${size.stripTrailingZeros().toPlainString()} ${units[unitIndex]}"
}

def formatPercent(value) {
    BigDecimal pct = parseDecimalMetric(value)
    if (pct == null) {
        return "NA"
    }
    return "${pct.setScale(1, java.math.RoundingMode.HALF_UP).stripTrailingZeros().toPlainString()}%"
}

def waitForTraceFile(File traceFile, long timeoutMs = 30000L) {
    long deadline = System.currentTimeMillis() + timeoutMs
    while (System.currentTimeMillis() <= deadline) {
        if (traceFile.exists() && traceFile.length() > 0L) {
            return true
        }
        sleep(200)
    }
    return traceFile.exists() && traceFile.length() > 0L
}

def readTraceRows(File traceFile) {
    if (!traceFile.exists() || traceFile.length() == 0L) {
        return []
    }

    def lines = traceFile.readLines().findAll { it?.trim() }
    if (lines.size() < 2) {
        return []
    }

    def headers = lines[0].split('\t', -1).collect { it.trim() }
    return lines.drop(1).collect { line ->
        def cols = line.split('\t', -1)
        def row = [:]
        headers.eachWithIndex { header, idx ->
            row[header] = idx < cols.size() ? cols[idx].trim() : ''
        }
        row.duration_ms = parseLongMetric(row['duration'])
        row.realtime_ms = parseLongMetric(row['realtime'])
        row.peak_rss_bytes = parseLongMetric(row['peak_rss'])
        row.peak_vmem_bytes = parseLongMetric(row['peak_vmem'])
        row.cpu_pct = parseDecimalMetric(row['%cpu'])
        row.cpus_req = parseDecimalMetric(row['cpus'])
        row.memory_req_bytes = parseLongMetric(row['memory'])
        row
    }
}

def padCell(value, int width) {
    def text = (value ?: 'NA').toString().replace('\t', ' ').replace('\n', ' ')
    if (text.size() > width) {
        if (width <= 3) {
            return text.take(width)
        }
        return text.take(width - 3) + '...'
    }
    return text.padRight(width)
}

def renderTable(List<Map> rows, List<List> columns) {
    if (!rows) {
        return ['(none)']
    }
    def header = columns.collect { col -> padCell(col[0], (int) col[2]) }.join('  ')
    def divider = columns.collect { col -> ''.padRight((int) col[2], '-') }.join('  ')
    def body = rows.collect { row ->
        columns.collect { col ->
            def key = col[1]
            def width = (int) col[2]
            padCell(row[key], width)
        }.join('  ')
    }
    return [header, divider] + body
}

def buildRunReport(workflowMeta, paramsMap) {
    def infoDir = new File(paramsMap.tracedir.toString())
    infoDir.mkdirs()

    def traceFile = new File(infoDir, 'trace.txt')
    def reportFile = new File(infoDir, 'run_stats.txt')
    def lines = []

    lines << 'BioPlumber Run Stats'
    lines << '===================='
    lines << ''
    lines << "roadmap_id: ${paramsMap.roadmap_id ?: 'NA'}"
    lines << "status: ${workflowMeta.success ? 'OK' : 'FAILED'}"
    lines << "run_name: ${workflowMeta.runName ?: 'NA'}"
    lines << "session_id: ${workflowMeta.sessionId ?: 'NA'}"
    lines << "started: ${workflowMeta.start ?: 'NA'}"
    lines << "completed: ${workflowMeta.complete ?: 'NA'}"
    lines << "elapsed: ${workflowMeta.duration ?: 'NA'}"
    lines << "command: ${workflowMeta.commandLine ?: 'NA'}"
    lines << "profile: ${workflowMeta.profile ?: 'NA'}"
    lines << "work_dir: ${workflowMeta.workDir ?: 'NA'}"
    lines << "trace_file: ${traceFile.absolutePath}"
    lines << ''

    def rows = waitForTraceFile(traceFile) ? readTraceRows(traceFile) : []
    if (!rows) {
        lines << 'Tasks'
        lines << '-----'
        lines << 'trace_status: unavailable or empty'
        lines << 'note: task-level CPU and memory stats depend on executor/container support'
        reportFile.text = lines.join(System.lineSeparator()) + System.lineSeparator()
        return
    }

    def totalTasks = rows.size()
    def completedTasks = rows.count { it.status == 'COMPLETED' }
    def failedTasks = rows.count { it.status in ['FAILED', 'ABORTED'] }
    def durations = rows.collect { it.duration_ms }.findAll { it != null }
    def peakRssValues = rows.collect { it.peak_rss_bytes }.findAll { it != null }
    def cpuValues = rows.collect { it.cpu_pct }.findAll { it != null }

    def processStats = rows.groupBy { it.process ?: 'unknown' }.collect { processName, processRows ->
        def processDurations = processRows.collect { it.duration_ms }.findAll { it != null }
        def processPeakRss = processRows.collect { it.peak_rss_bytes }.findAll { it != null }
        def processCpu = processRows.collect { it.cpu_pct }.findAll { it != null }
        [
            process: processName,
            tasks: processRows.size(),
            total_duration: processDurations ? processDurations.inject(0L) { acc, v -> acc + v } : 0L,
            max_duration: processDurations ? processDurations.max() : null,
            max_peak_rss: processPeakRss ? processPeakRss.max() : null,
            max_cpu: processCpu ? processCpu.max() : null
        ]
    }.sort { a, b ->
        def durationCmp = (b.total_duration ?: 0L) <=> (a.total_duration ?: 0L)
        durationCmp != 0 ? durationCmp : (a.process <=> b.process)
    }

    def slowestProcess = processStats ? processStats[0].process : 'NA'
    def highestMemoryProcess = processStats.findAll { it.max_peak_rss != null }.sort { a, b ->
        (b.max_peak_rss ?: 0L) <=> (a.max_peak_rss ?: 0L)
    }
    highestMemoryProcess = highestMemoryProcess ? highestMemoryProcess[0].process : 'NA'

    lines << 'Tasks'
    lines << '-----'
    lines << "total_tasks: ${totalTasks}"
    lines << "completed_tasks: ${completedTasks}"
    lines << "failed_tasks: ${failedTasks}"
    lines << ''
    lines << 'Peaks'
    lines << '-----'
    lines << "max_task_duration: ${durations ? formatDurationMs(durations.max()) : 'NA'}"
    lines << "max_peak_rss: ${peakRssValues ? formatBytes(peakRssValues.max()) : 'NA'}"
    lines << "max_percent_cpu: ${cpuValues ? formatPercent(cpuValues.max()) : 'NA'}"
    lines << "slowest_process_by_total_duration: ${slowestProcess}"
    lines << "highest_memory_process: ${highestMemoryProcess}"
    lines << ''
    lines << 'Top Processes By Total Duration'
    lines << '-------------------------------'

    def topProcesses = processStats.take(10).collect { row ->
        [
            process: row.process,
            tasks: row.tasks,
            total_duration: formatDurationMs(row.total_duration),
            max_duration: formatDurationMs(row.max_duration),
            max_peak_rss: formatBytes(row.max_peak_rss),
            max_cpu: formatPercent(row.max_cpu)
        ]
    }

    lines.addAll(renderTable(
        topProcesses,
        [
            ['process', 'process', 32],
            ['tasks', 'tasks', 7],
            ['total_duration', 'total_duration', 16],
            ['max_duration', 'max_duration', 16],
            ['max_peak_rss', 'max_peak_rss', 14],
            ['max_%cpu', 'max_cpu', 10]
        ]
    ))

    lines << ''
    lines << 'Top Tasks By Duration'
    lines << '---------------------'

    def topTasks = rows.findAll { it.duration_ms != null }.sort { a, b ->
        (b.duration_ms ?: 0L) <=> (a.duration_ms ?: 0L)
    }.take(10).collect { row ->
        [
            process: row.process ?: 'NA',
            task: row.tag ?: row.name ?: 'NA',
            duration: formatDurationMs(row.duration_ms),
            peak_rss: formatBytes(row.peak_rss_bytes),
            cpu: formatPercent(row.cpu_pct),
            status: row.status ?: 'NA'
        ]
    }

    lines.addAll(renderTable(
        topTasks,
        [
            ['process', 'process', 28],
            ['task', 'task', 28],
            ['duration', 'duration', 14],
            ['peak_rss', 'peak_rss', 12],
            ['%cpu', 'cpu', 9],
            ['status', 'status', 10]
        ]
    ))

    lines << ''
    lines << 'Notes'
    lines << '-----'
    lines << 'trace.txt contains one row per task execution.'
    lines << 'The task column above prefers the Nextflow tag when available, then falls back to the task name.'
    lines << 'CPU and memory metrics depend on executor/container support; missing values are reported as NA.'

    reportFile.text = lines.join(System.lineSeparator()) + System.lineSeparator()
}

workflow {
    if (params.roadmap_id=="roadmap_1")
    {

            if (params.input_type=="sra")
            {
                table=tableToDict(file("${params.input_file}"))
                get_sequences_from_sra(Channel.fromList(table["Run"]))
                roadmap_1(get_sequences_from_sra.out.sra_ids, get_sequences_from_sra.out.fastq_files)
            }
            if (params.input_type=="local")
            {
                table=tableToDict(file("${params.input_file}"))
                if (!params.single_end)
                {
                reads_1=Channel.fromPath(table["reads1"])
                reads_2=Channel.fromPath(table["reads2"])
                reads=reads_1.merge(reads_2)  
                }
                else
                {
                reads=Channel.fromPath(table["reads"].collect{t->file(t)})
                }

                sample_name=Channel.fromList(table["sample_name"])
                roadmap_1(sample_name, reads)
            }
            }
        
        else if (params.roadmap_id=="roadmap_2")
        {
            if (params.input_reads)
            {
                table=tableToDict(file("${params.input_reads}"))
                reads_1=Channel.fromPath(table["reads1"].collect{t->file(t)})
                reads_2=Channel.fromPath(table["reads2"].collect{t->file(t)})
                inputs=reads_1.merge(reads_2)
                sample_names=Channel.from(table["sample_name"])
                params.roadmap_2_input_type="reads"
            }
            else if (params.input_bams)
            {
                table=tableToDict(file("${params.input_bams}"))
                inputs=Channel.fromPath(table["bam_files"].collect{t->file(t)})
                sample_names=Channel.from(table["sample_name"])
                params.roadmap_2_input_type="bams"
            }
            else
            {
                error "Please provide the reads or bam files information using the input_reads or input_bam parameter."
            }
            
            
            
            if (params.input_fastas)
            {
            if (params.is_genome_db||params.is_stb_db)
            {
                error "If you provide a list of fasta files, you cannot provide a genome database or STB file."
            }
            else
            {
            fasta_file=tableToDict(file("${params.input_fastas}"))["fasta_files"].collect{t->file(t)}
            }
            }

            if (params.is_genome_db)
            {
                if (!params.is_stb_db)
                {
                    error "If you provide a genome database, you must also provide an STB file."
                }
                else{
                    fasta_file=file(params.is_genome_db)
                }
            }
            roadmap_2(sample_names, inputs, fasta_file)
        }



        else if (params.roadmap_id=="roadmap_3")
        {
            if (params.input_type=="path")
            {
                genomes=files("${params.input_genomes}")

                roadmap_3(genomes)
            }
        }
        else if (params.roadmap_id=="roadmap_1_3_2")
        {
            if (!params.host_genome)
            {
                error "Please provide a host genome for decontamination in the roadmap_1 workflow."
            }
            params.roadmap_2_input_type="reads"
            if (params.input_type=="sra")
            {
                table=tableToDict(file("${params.input_file}"))
                get_sequences_from_sra(Channel.fromList(table["Run"]))
                sample_names=get_sequences_from_sra.out.sra_ids
                reads=get_sequences_from_sra.out.fastq_files
                roadmap_1_3_2(sample_names, reads, file(params.host_genome))
            }
            if (params.input_type=="local")
            {
                table=tableToDict(file("${params.input_file}"))
                reads_1=Channel.fromPath(table["reads1"])
                reads_2=Channel.fromPath(table["reads2"])
                reads=reads_1.merge(reads_2)
                sample_name=table["sample_name"]
                sample_name=Channel.fromList(sample_name)
                roadmap_1_3_2(sample_name, reads, file(params.host_genome))
            }
        }
    else if (params.roadmap_id=="roadmap_4")
    {
        if (!params.host_genome)
            {
                error "Please provide a host genome for decontamination in the roadmap_1 workflow."
            }


            if (params.input_type=="sra")
            {
                table=tableToDict(file("${params.input_file}"))
                get_sequences_from_sra(Channel.fromList(table["Run"]))
                roadmap_4(get_sequences_from_sra.out.sra_ids, get_sequences_from_sra.out.fastq_files, file(params.host_genome))
            }
            if (params.input_type=="local")
            {
                table=tableToDict(file("${params.input_file}"))
                reads_1=Channel.fromPath(table["reads1"])
                reads_2=Channel.fromPath(table["reads2"])
                reads=reads_1.merge(reads_2)
                sample_name=Channel.fromList(table["sample_name"])
                roadmap_4(sample_name, reads, file(params.host_genome))
            }
    }

    else if (params.roadmap_id=="roadmap_3_2")
    {
            if (!params.input_reads)
            {
                error "Please provide the reads information using the input_reads parameter."
            }
            if (!params.input_fastas)
            {
                error "Please provide the fasta files information using the input_fastas parameter."
            }
            table=tableToDict(file("${params.input_reads}"))
            reads_1=Channel.fromPath(table["reads1"].collect{t->file(t)})
            reads_2=Channel.fromPath(table["reads2"].collect{t->file(t)})
            reads=reads_1.merge(reads_2)
            sample_names=Channel.from(table["sample_name"])
            fasta_file=tableToDict(file("${params.input_fastas}"))["fasta_files"].collect{t->file(t)}
            if (params.force_genomes)
            {
              force_genomes=Channel.fromPath(file(params.force_genomes))  
            }
            else
            {
              force_genomes=Channel.empty()
            }
            params.roadmap2_input_type="reads"

            roadmap_3_2(sample_names, reads, fasta_file, force_genomes)
    
    }
    else if (params.roadmap_id=="roadmap_5")
    {
        if (!params.input_fastas)
        {
            error "Please provide a genome for mapping reads."
        }
        if (!params.input_reads)
        {
            error "Please provide the reads information using the input_reads parameter."
        }
        def valid_read_types = ["short", "nanopore", "pacbio_clr", "pacbio_hifi"]
        if (!valid_read_types.contains(params.read_type))
        {
            error "Invalid --read_type '${params.read_type}'. Must be one of: ${valid_read_types.join(', ')}"
        }
        if (params.read_type == "short")
        {
            // Short reads: local CSV needs columns: sample_name, reads1, reads2; SRA CSV needs column: Run
            if (params.input_type=="sra")
            {
                table=tableToDict(file("${params.input_reads}"))
                get_sequences_from_sra(Channel.fromList(table["Run"]))
                sample_names=get_sequences_from_sra.out.sra_ids
                reads=get_sequences_from_sra.out.fastq_files
            }
            else if (params.input_type=="local")
            {
                table=tableToDict(file("${params.input_reads}"))
                reads_1=Channel.fromPath(table["reads1"].collect{t->file(t)})
                reads_2=Channel.fromPath(table["reads2"].collect{t->file(t)})
                reads=reads_1.merge(reads_2)
                sample_names=Channel.fromList(table["sample_name"])
            }
            samples_reads=sample_names.merge(reads)
            genomes=Channel.fromPath(tableToDict(file("${params.input_fastas}"))["fasta_files"].collect{t->file(t)})
            if (params.roadmap_5_pairmode=="cross")
            {
                samples_reads.combine(genomes).set{inputs}
            }
            else if (params.roadmap_5_pairmode=="paired")
            {
                samples_reads.merge(genomes).set{inputs}
            }
            inputs.multiMap{v->
                sn:v[0]
                rd:v[1]
                gn:v[2]
                tr:(v[1].size() == 2)
            }.set{ins}
            roadmap_5(ins.sn, ins.rd, ins.gn, ins.tr)
        }
        else
        {
            // Long reads (nanopore, pacbio_clr, pacbio_hifi): CSV must have columns: sample_name, reads (local) or Run (sra)
            if (params.input_type=="sra")
            {
                table=tableToDict(file("${params.input_reads}"))
                get_sequences_from_sra(Channel.fromList(table["Run"]))
                sample_names=get_sequences_from_sra.out.sra_ids
                reads=get_sequences_from_sra.out.fastq_files
            }
            else if (params.input_type=="local")
            {
                table=tableToDict(file("${params.input_reads}"))
                reads=Channel.fromPath(table["reads"].collect{t->file(t)})
                sample_names=Channel.fromList(table["sample_name"])
            }
            samples_reads=sample_names.merge(reads)
            genomes=Channel.fromPath(tableToDict(file("${params.input_fastas}"))["fasta_files"].collect{t->file(t)})
            if (params.roadmap_5_pairmode=="cross")
            {
                samples_reads.combine(genomes).set{inputs}
            }
            else if (params.roadmap_5_pairmode=="paired")
            {
                samples_reads.merge(genomes).set{inputs}
            }
            inputs.multiMap{v->
                sn: v[0]
                rd: v[1]
                gn: v[2]
                tr: false   // long reads are always single-end
            }.set{ins}
            roadmap_5(ins.sn, ins.rd, ins.gn, ins.tr)
        }
    }
    else if (params.roadmap_id=="roadmap_6")
    {
        if (params.input_type=="sra")
        {
            table=tableToDict(file("${params.input_file}"))
            get_sequences_from_sra(Channel.fromList(table["Run"]))
            sample_names=get_sequences_from_sra.out.sra_ids
            reads=get_sequences_from_sra.out.fastq_files
            roadmap_6(sample_names, reads)
        }
        else if (params.input_type=="local")
        {
            table=tableToDict(file("${params.input_file}"))
            reads_1=Channel.fromPath(table["reads1"].collect{t->file(t)})
            reads_2=Channel.fromPath(table["reads2"].collect{t->file(t)})
            reads=reads_1.merge(reads_2)
            sample_name=Channel.fromList(table["sample_name"])
            roadmap_6(sample_name, reads)
        }
        else
        {
            error "Please provide the reads information using the file parameter."
        }
    }
    else if (params.roadmap_id=="roadmap_7")
    {
        if (params.bins_dir)
        {
            bins=Channel.fromPath(params.bins_dir,type: 'file')
        }
        else if(params.input_bins_table)
        {
            table=tableToDict(file("${params.input_bins_table}"))
            bins=Channel.fromPath(table["fasta_files"].collect{t->file(t)})

        }
        else
        {
            error "Please provide the bins information using either the bins_dir or input_bins_table parameter."
        }
        roadmap_7(bins.collect())
        
    }
    else if (params.roadmap_id=="roadmap_8")
    {
        if (params.input_type=="sra")
        {
            table=tableToDict(file("${params.input_file}"))
            get_sequences_from_sra(Channel.fromList(table["Run"]))
            sample_names=get_sequences_from_sra.out.sra_ids
            reads=get_sequences_from_sra.out.fastq_files
            host_genome=file(params.host_genome)
            host_genome_gtf=file(params.host_genome_gtf)
        }
        else if (params.input_type=="local")
        {
            table=tableToDict(file("${params.input_file}"))
            reads_1=Channel.fromPath(table["reads1"].collect{t->file(t)})
            reads_2=Channel.fromPath(table["reads2"].collect{t->file(t)})
            reads=reads_1.merge(reads_2)
            sample_name=Channel.fromList(table["sample_name"])
            host_genome=file(params.host_genome)
            host_genome_gtf=file(params.host_genome_gtf)

        }

        else
        {
            error "Please provide the reads information using the file parameter."
        }

        if (params.mode=="bulk_rna_seq")
        {
            bulk_rna_seq(sample_names, reads, host_genome, host_genome_gtf)
        }
        else if (params.mode=="single_cell_rna_seq")
        {   
            single_cell_rna_seq(sample_names, reads, host_genome, host_genome_gtf)
        }
        else
        {
            error "Please provide a valid mode."
        }
    }
    else if (params.roadmap_id=="roadmap_9")
    {
        if (params.input_type=="sra")
        {
            table=tableToDict(file("${params.input_file}"))
            get_sequences_from_sra(Channel.fromList(table["Run"]))
            sample_names=get_sequences_from_sra.out.sra_ids
            reads=get_sequences_from_sra.out.fastq_files
            reference_transcriptome=file(params.reference_transcriptome)
        }
        else if (params.input_type=="local")
        {
            table=tableToDict(file("${params.input_file}"))
            reads_1=Channel.fromPath(table["reads1"].collect{t->file(t)})
            if (params.paired_end)
            {
            reads_2=Channel.fromPath(table["reads2"].collect{t->file(t)})
            reads=reads_1.merge(reads_2)
            }
            else{
                reads_1.set{reads}
            }

            sample_names=Channel.fromList(table["sample_name"])
            reference_genome=file(params.reference_genome)

        }
        else
        {
            error "Please provide the reads information using the file parameter."
        }
        roadmap_9(sample_names, reads, reference_genome)

    }
    else if (params.roadmap_id=="quality_control")
    {
        if (!params.host_genome)
            {
                error "Please provide a host genome for decontamination in the QC workflow."
            }
            if (params.input_type=="sra")
            {
                table=tableToDict(file("${params.input_file}"))
                get_sequences_from_sra(Channel.fromList(table["Run"]))
                quality_control(get_sequences_from_sra.out.sra_ids, get_sequences_from_sra.out.fastq_files, file(params.host_genome))
            }
            if (params.input_type=="local")
            {
                table=tableToDict(file("${params.input_file}"))
                reads_1=Channel.fromPath(table["reads1"])
                reads_2=Channel.fromPath(table["reads2"])
                reads=reads_1.merge(reads_2)
                sample_name=Channel.fromList(table["sample_name"])
                quality_control(sample_name, reads, file(params.host_genome))
            }
    }
    else if (params.roadmap_id=="annotate_contigs")
    {
        if (!params.input_contigs)
        {
            error "Please provide the contigs information using the input_contigs parameter."
        }
        table=tableToDict(file("${params.input_contigs}"))
        sample_names=Channel.fromList(table["sample_name"])
        contigs=Channel.fromList(table["contig_files"].collect{t->file(t)})
        annotate_contigs(sample_names, contigs)
    }

    else if (params.roadmap_id=="download_samples")
    {

        table=tableToDict(file("${params.input_file}"))
        get_sequences_from_sra(Channel.fromList(table["Run"]))

    }
    else if (params.roadmap_id=="download_links")
    {

        table=tableToDict(file("${params.input_file}"))
        links=Channel.fromList(table["Links"])
        names=Channel.fromList(table["Names"])
        all=links.merge(names)
        download_files(all)

    }
    else if (params.roadmap_id=="scan_metagenome")
    {

        table=tableToDict(file("${params.input_file}"))
        Channel.fromList(table["Run"]).set{ sra_acc }
        genome=file(params.genome)
        scanGenome(sra_acc.collate(params.scan_genome_batch_size), genome)


    }
    
    else if (params.roadmap_id=="subsample_reads")
    {
        if (!params.input_file)
        {
            error "Please provide the reads information using the --input_file parameter."
        }
        if (!params.fractions)
        {
            error "Please provide fractions using the --fractions parameter (e.g. --fractions '0.1,0.25,0.5')."
        }
        fractions=Channel.fromList(params.fractions.tokenize(',').collect{it.trim() as Double})
        if (params.input_type=="sra")
        {
            table=tableToDict(file("${params.input_file}"))
            get_sequences_from_sra(Channel.fromList(table["Run"]))
            sample_names=get_sequences_from_sra.out.sra_ids
            reads=get_sequences_from_sra.out.fastq_files
        }
        else if (params.input_type=="local")
        {
            table=tableToDict(file("${params.input_file}"))
            reads_1=Channel.fromPath(table["reads1"].collect{t->file(t)})
            reads_2=Channel.fromPath(table["reads2"].collect{t->file(t)})
            reads=reads_1.merge(reads_2)
            sample_names=Channel.fromList(table["sample_name"])
        }
        else
        {
            error "Please provide a valid --input_type (local or sra)."
        }
        subsample_reads(sample_names, reads, fractions)
    }
    else
        {
            error "Please provide a valid roadmap_id."
        }
        
        
    }

    

// ###### Roadmaps ###### //


workflow roadmap_1{

    take:
    sample_name
    reads

    main:
    
    assembly(sample_name, reads)
    binning(assembly.out.sample_name, assembly.out.sorted_bams, assembly.out.contigs,assembly.out.reads)
    emit:
    metabat2_bins=binning.out.metabat2_bins
}


workflow roadmap_2 {
    // This roadmap takes a list of reads, their corresponding sample names, and a list of fasta files.
    // It performs strain-level profiling using InStrain.

    take:
    sample_names
    inputs
    fasta_file
    main:
    if (!params.is_genome_db)
    {
        if (params.add_fasta_prefix)
        {
            add_prefix_to_fasta(fasta_file)
            fasta_file=add_prefix_to_fasta.out.prefixed_fasta
        }
        else
        {
            fasta_file=fasta_file
        }
        concatenate_files(fasta_file, "genomes_db.fasta")
        genome_db=concatenate_files.out.concatenated_file
    }
    else
    {
        genome_db=file(params.is_genome_db)
    }

    if (!params.is_stb_db)
    {

        make_stb_file_instrain(fasta_file, "genomes_db")
        stb=make_stb_file_instrain.out.stb_file
    }
    else
    {
        stb=file(params.is_stb_db)
    }

    if (!params.is_genes)
    {
        find_genes_prodigal(genome_db)
        genes=find_genes_prodigal.out.genes_fna
    }
    else
    {
        genes=file(params.is_genes)
    }

    if (params.roadmap_2_input_type=="reads")
    {
        index_bowtie2(genome_db, "genomes_db")
        bowtie2_to_sorted_bam(sample_names, index_bowtie2.out.reference_genome, inputs, index_bowtie2.out.bowtie2_index_files,false)
        sorted_bam=bowtie2_to_sorted_bam.out.sorted_bam
    }
    else if (params.roadmap_2_input_type=="bams")
    {
        sorted_bam=inputs
    }
    profile_with_instrain(sorted_bam,genome_db, stb, genes)
    profile_with_instrain.out.instrain_profiles.collect().set{all_profiles}
    compare_instrain_profiles(all_profiles, stb)
    emit:
    instrain_profiles=all_profiles


}
    
workflow roadmap_3 {
    take:
    genomes
    
    main:
    write_genome_list(genomes)
    dereplicate_drep(genomes,write_genome_list.out.genomes_list)

    emit:
    dereplicated_genomes=dereplicate_drep.out.dereplicated_genomes

}

workflow roadmap_4 {
    take:
    sample_name
    reads
    host_genome
    main:
    quality_control(sample_name, reads, host_genome)
    emit:
    qc_reads=quality_control.out.qc_reads
    
}

workflow roadmap_1_3_2{
    /*
    This workflow gets the bins from individual samples using roadmap 1, dereplicates them using roadmap 3, and then performs strain-level profiling using roadmap 2.
    WARNING: This workflow is not fully tested and may not work as expected.
    */
    take:
    sample_name
    reads
    host_genome

    main:
    roadmap_1(sample_name, reads, host_genome)
    params.genomes_exctension=params.binning_extension
    dereplicated_genomes=roadmap_3(roadmap_1.out.metabat2_bins.collect())
    params.is_genome_db=null
    params.is_stb_db=null
    roadmap_2(sample_name, reads, dereplicated_genomes)


}

workflow roadmap_3_2 {
    // This roadmap is appropriate for ddereplicating a series of bins and then performing strain-level profiling using InStrain on the dereplicated bins.
    // An example would be comparative genomics for closely related bins
    take:
    sample_names
    reads
    fasta_file
    external_genomes
    main:
    roadmap_3(fasta_file.collect())
    dereplicated_genomes=roadmap_3.out.dereplicated_genomes
    dereplicate_drep=dereplicated_genomes.mix(external_genomes).flatten().unique().collect()
    roadmap_2(sample_names, reads, dereplicate_drep)

    emit:
    instrain_profiles=roadmap_2.out.instrain_profiles

}

workflow roadmap_5 {
    /*
     * Maps reads to reference genomes. Supports short reads (bowtie2) and long reads
     * (minimap2: nanopore, pacbio_clr, pacbio_hifi). Controlled by --read_type.
     *
     * Input modes (--input_type):
     *   sra   : CSV with a "Run" column containing SRA accession IDs. Supported for
     *           both short and long reads.
     *   local : CSV with "sample_name" + "reads1"/"reads2" columns (short reads) or
     *           "sample_name" + "reads" columns (long reads).
     *
     * Required parameters:
     *   --input_reads  : path to the input CSV file (reads or SRA accessions)
     *   --input_fastas : path to a CSV with a "fasta_files" column listing reference genomes
     *   --read_type    : short | nanopore | pacbio_clr | pacbio_hifi (default: short)
     *
     * Pairing mode (--roadmap_5_pairmode):
     *   paired : each sample is mapped to the genome at the same row in input_fastas (default)
     *   cross  : every sample is mapped to every genome (all-vs-all)
     *
     * Output options:
     *   --keep_unmapped_reads : output BAM retains unmapped reads as well (default: mapped-only)
     *   --get_mapped_reads    : additionally produce a FASTQ of mapped reads
     *   --get_unmapped_reads  : additionally produce a FASTQ of unmapped reads
     */
    take:
    sample_name
    reads
    genome
    paired

    main:
    if (params.read_type == "short") {
        map_reads_fasta_pairs(sample_name, reads, genome, paired)
    } else {
        map_long_reads_fasta_pairs(sample_name, reads, genome)
    }
}

workflow subsample_reads {
    /*
     * Randomly subsamples reads at one or more fractions using BBTools reformat.sh.
     * Paired-end reads are processed together so read pairing is preserved.
     * Each sample × fraction combination produces its own output FASTQ file(s).
     *
     * Input modes (--input_type):
     *   sra   : CSV with a "Run" column containing SRA accession IDs.
     *   local : CSV with "sample_name", "reads1", "reads2" columns.
     *
     * Required parameters:
     *   --input_file : path to the input CSV file
     *   --fractions  : comma-separated list of fractions (e.g. "0.1,0.25,0.5")
     *
     * Optional parameters:
     *   --subsample_seed : random seed for reproducibility (default: 42)
     */
    take:
    sample_names
    reads
    fractions

    main:
    sample_names.merge(reads)
        .combine(fractions)
        .multiMap { v ->
            sn: v[0]
            rd: v[1]
            fr: v[2]
        }.set { ins }
    subsample_reads_reformat(ins.sn, ins.rd, ins.fr)
}

workflow roadmap_6{
    // This roadmap analyzes QCed metagenomics reads and performs functional and taxonomic profiling using
    // Reference-based methods.
    take:
    sample_name
    reads
    main:

    sample_name.multiMap{t->
        sample_name_sylph:t
        sample_name_metaphlan:t
    }.set{sample_names}


    reads.multiMap{t->
        reads_sylph:t
        reads_metaphlan:t
    }.set{reads_all}

    reads_all.reads_sylph.multiMap{t->
    reads1_sylph: t[0]
    reads2_sylph: t[1]
    }.set{sylph_reads}

    if (!params.exclude_sylph)
    {
        if (params.sylph_db)
    {
        sylph_db=file(params.sylph_db)
    }
    else
    {
        sylph_db=download_sylph_db()
    }

    estimate_abundance_sylph(sylph_reads.reads1_sylph.collect(), sylph_reads.reads2_sylph.collect(), sylph_db)
    }
   
    if (!params.exclude_metaphlan)
    {
        if (params.metaphlan_db)
          {
           metaphlan_db=file(params.metaphlan_db)
          }
        else
          {
          metaphlan_db=download_metaphlan_db()
          }
        estimate_abundance_metaphlan(sample_names.sample_name_metaphlan,reads_all.reads_metaphlan, metaphlan_db)
        merge_metaphlan_tables(estimate_abundance_metaphlan.out.abundance.collect())
        calculate_diversity_metaphlan(merge_metaphlan_tables.out.merged_table)
    }
   
    if (!params.exclude_kraken)
    {
        if (params.kraken2_db)
        {
            kraken2_db=file(params.kraken2_db)
        }
        else
        {
            kraken2_db=download_kraken2_db(params.kraken2_db_link)
        }
        classify_kraken2(sample_name,reads, kraken2_db)
        estimate_abundance_bracken(classify_kraken2.out.sample_name,classify_kraken2.out.kraken_report,kraken2_db)
    }

    
    if (!params.exclude_humann)
    {
        if (params.humann_uniref90)
        {
            humann_uniref90=file(params.humann_uniref90)
        }
        else
        {
	    download_humann_uniref90()
            humann_uniref90=download_humann_uniref90.out.humann_uniref90
        }
        if (params.humann_chocophlan)
        {
            humann_chocophlan=file(params.humann_chocophlan)
        }
        else
        {
            download_humann_chocophlan()
	    humann_chocophlan=download_humann_chocophlan.out.humann_chocophlan
        }
        profile_humann(sample_name, reads, humann_chocophlan, humann_uniref90, metaphlan_db)    
    }
    if (!params.exclude_eukdetect)
    {
        if (params.eukdetect_db)
        {
            eukdetect_db=file(params.eukdetect_db)
        }
        else
        {
            error "Please provide a database for EukDetect. This is required for EukDetect to run."
        }
        eukdetect(sample_name, reads, eukdetect_db)
    }
}

workflow roadmap_7{
    // This roadmap takes a set of bins and performs taxonomic assignment using GTDB-Tk.//
    take:
    bins
    main:
    if (params.gtdbtk_db)
    {
        gtdbtk_db=file(params.gtdbtk_db)
    }
    else
    {
        gtdbtk_db=download_gtdbtk_db()
    }
    assign_taxonomy_gtdb_tk(bins, gtdbtk_db)
    

}

workflow roadmap_9{
    take:
    sample_names
    reads
    reference_genome
    main:
        assemble_rna_spades(sample_names, reads)
        contigs=assemble_rna_spades.out.soft_filtered_transcripts.concat(assemble_rna_spades.out.hard_filtered_transcripts)
        sample_names_ordered=assemble_rna_spades.out.sample_name.concat(assemble_rna_spades.out.sample_name)
        contigs.merge(sample_names_ordered).set{contigs_with_sn}
        contigs_with_sn.multiMap{t->
            sample_name: t[1]+"_"+t[0].baseName
            contig_file:t[0]
            }.set{contigs_to_be_processed}
        get_circular_contigs_cirit(contigs_to_be_processed.sample_name, contigs_to_be_processed.contig_file)
        map_rna_assemblies_to_reference_genome(get_circular_contigs_cirit.out.circular_contigs.collect(),reference_genome)

    emit:
    mapped_contigs=map_rna_assemblies_to_reference_genome.out.mapped_contigs
    unmapped_contigs=map_rna_assemblies_to_reference_genome.out.unmapped_contigs

}

workflow bulk_rna_seq{
    // Standard bulk RNA-Seq workflow for quality control, alignment, and quantification.

    take:
    sample_name
    reads
    host_genome
    host_genome_gtf
    main:
    read_qc_fastp(sample_name, reads)
    index_star(host_genome, host_genome_gtf)
    align_star(read_qc_fastp.out.sample_name,index_star.out.star_index_files, read_qc_fastp.out.fastp_qcd_reads)
    align_star.out.star_aligned_bam.collect().set{all_bams}
    gene_count_featurecounts(all_bams, host_genome_gtf)


}
// ###### WORKFLOWS ###### //
workflow single_cell_rna_seq{
    take:
    sample_name
    reads
    host_genome
    host_genome_gtf 
    main:
    read_qc_fastp(sample_name, reads)
    index_kallisto(host_genome, host_genome_gtf)
    map_reads_kallisto_single_cell(sample_name, index_kallisto.out.index_file, index_kallisto.out.t2g_file, read_qc_fastp.out.fastp_qcd_reads)
}
workflow quality_control {
    /*
    This workflow takes raw sequencing reads and performs quality control and decontamination.
    It generates cleaned reads and quality reports.
    */
    take:
    sample_name
    reads
    host_genome

    main:
    read_qc_fastp(sample_name, reads)
    index_bowtie2(host_genome,host_genome.baseName)
    bowtie2_to_sorted_bam(read_qc_fastp.out.sample_name,index_bowtie2.out.reference_genome, read_qc_fastp.out.fastp_qcd_reads,index_bowtie2.out.bowtie2_index_files,true)  
    get_unmapped_reads(bowtie2_to_sorted_bam.out.sorted_bam,bowtie2_to_sorted_bam.out.paired,bowtie2_to_sorted_bam.out.sample_name)
    get_mapped_reads(bowtie2_to_sorted_bam.out.sorted_bam,bowtie2_to_sorted_bam.out.paired,bowtie2_to_sorted_bam.out.sample_name)

    emit:
        qc_reads=get_unmapped_reads.out.unmapped_reads
        sample_name=get_unmapped_reads.out.sample_name
        paired=bowtie2_to_sorted_bam.out.paired
}

workflow scanGenome {
    take:
    sra_acc
    genome
    main:
    buildSylphdb(genome)
    sylphScanBatchFromSRA(sra_acc, buildSylphdb.out.sylph_db)
}

workflow assembly {
    take:
    sample_name
    reads

    main:
    assemble_with_megahit(sample_name, reads)

    map_reads_fasta_pairs(   
                             assemble_with_megahit.out.sample_name,
                             assemble_with_megahit.out.reads,
                             assemble_with_megahit.out.contigs,
                             assemble_with_megahit.out.paired
                        )
    
    emit:
    contigs=map_reads_fasta_pairs.out.reference_fasta
    sorted_bams=map_reads_fasta_pairs.out.sorted_bam
    sample_name=map_reads_fasta_pairs.out.sample_name
    reads=map_reads_fasta_pairs.out.reads
}
    
workflow binning {
    take:
    sample_name
    sorted_bams
    assembly
    reads
    main:
    get_coverage_for_metabat2(sample_name, sorted_bams, assembly)
    binning_with_metabat2(get_coverage_for_metabat2.out.sample_name,get_coverage_for_metabat2.out.contig,get_coverage_for_metabat2.out.coverage)

    emit:
    metabat2_bins=binning_with_metabat2.out.metabat2_bins
    sample_name=binning_with_metabat2.out.sample_name
    reads=reads

}

workflow test_customized_compared{
    take:
    profile_pairs
    stb_file
    main:
    compare_general_customized(profile_pairs,stb_file)

}


workflow annotate_contigs{
    take:
    sample_names
    contigs
    main:
    find_genes_prodigal(contigs)
    if (params.kraken2_db)
        {
            kraken2_db=file(params.kraken2_db)
        }
    else
        {
            kraken2_db=download_kraken2_db(params.kraken2_db_link)
        }
    classify_kraken2_contigs(sample_names, contigs, kraken2_db)

    if (params.build_gene_db_mode=="nucleotide")
        {
        concatenate_files(find_genes_prodigal.out.genes_fna.collect(), "all_genes.fna")
        }
    else if (params.build_gene_db_mode=="amino_acid")
        {
        concatenate_files(find_genes_prodigal.out.genes_faa.collect(), "all_genes.faa")
        }
    else
    {
        error "Please provide a valid build_gene_db_mode parameter: nucleotide or amino_acid."
    }
    if (!params.skip_mmseqs_clustering)
        {
        mmseqs_linclust(concatenate_files.out.concatenated_file)
        }
    if (!params.skip_functional_annotation)
    {
        if (params.eggnog_data_dir)
        {
            egg_db=file(params.eggnog_data_dir)
        }
        else
        {
            download_eggnog_db()
            egg_db=download_eggnog_db.out.eggnog_db
        }
        if (!params.skip_mmseqs_clustering)
        {
            eggnog_annotation(mmseqs_linclust.out.clustered_sequences, egg_db)
        }
        else
        {
            eggnog_annotation(concatenate_files.out.concatenated_file, egg_db)
        }
    }
    if (!params.skip_genomad_annotation)
    {
        if (params.genomad_db)
        {
            genomad_db=file(params.genomad_db)
        }
        else
        {
            download_genomad_db()
            genomad_db=download_genomad_db.out.genomad_db
        }
        annotate_contig_genomad(sample_names, contigs, genomad_db)
    }
    
}





// Include modules

include {read_qc_fastp} from "./modules/qc"

include {
    index_bowtie2;
    align_bowtie2;
    convert_sam_to_bam;
    sort_bam;
    get_unmapped_reads;
    get_mapped_reads;
    map_reads_fasta_pairs;
    map_long_reads_fasta_pairs;
    bowtie2_to_sorted_bam;
    index_star;
    align_star;
    index_kallisto;
    map_reads_kallisto_single_cell;
    map_contigs_to_reference_transcriptome;
    map_rna_assemblies_to_reference_genome;
        } from "./modules/alignment"

include {assemble_with_megahit;
        assemble_rna_spades;
        get_circular_contigs_cirit;
        } from "./modules/assembly"

include {get_coverage_for_metabat2;
        binning_with_metabat2;
        add_prefix_to_fasta} from "./modules/binning"

include {tableToDict;
        concatenate_files;
        subsample_reads_reformat;
        } from "./modules/files"

include {get_sequences_from_sra;
         download_files;
} from "./modules/download"

include {estimate_abundance_coverm;
         estimate_abundance_metaphlan;
            estimate_abundance_sylph;
            merge_metaphlan_tables;
            calculate_diversity_metaphlan;
            download_sylph_db;
            download_metaphlan_db;
            download_kraken2_db;
            classify_kraken2;
            estimate_abundance_bracken;
            gene_count_featurecounts;
            sylphScanBatchFromSRA;
            buildSylphdb;
            eukdetect;

         } from './modules/abundance'

include {compare_instrain_profiles;
         profile_with_instrain;
         make_stb_file_instrain;
         compare_general_customized;
         sample_pairs;
         } from './modules/strain'

include { dereplicate_drep;write_genome_list } from './modules/dereplication'

include {find_genes_prodigal} from './modules/genes'

include {
    download_humann_chocophlan;
    download_humann_uniref90;
    profile_humann;
} from './modules/metabolism'

include {assign_taxonomy_gtdb_tk;
        download_gtdbtk_db;
        classify_kraken2_contigs;
        eggnog_annotation;
        download_eggnog_db;
        download_genomad_db;
        annotate_contig_genomad;
        } from './modules/annotation'

include {
 create_mmseqs_db;
 mmseqs_linclust;
} from './modules/cluster'

workflow.onComplete {
    if (!params.dump_report) {
        return
    }

    try {
        buildRunReport(workflow, params)
    } catch (Exception e) {
        def infoDir = new File(params.tracedir.toString())
        infoDir.mkdirs()
        def reportFile = new File(infoDir, 'run_stats.txt')
        reportFile.text = """BioPlumber Run Stats
====================

roadmap_id: ${params.roadmap_id ?: 'NA'}
status: ${workflow.success ? 'OK' : 'FAILED'}
run_name: ${workflow.runName ?: 'NA'}
session_id: ${workflow.sessionId ?: 'NA'}
started: ${workflow.start ?: 'NA'}
completed: ${workflow.complete ?: 'NA'}
elapsed: ${workflow.duration ?: 'NA'}
trace_file: ${new File(params.tracedir.toString(), 'trace.txt').absolutePath}
report_generation_error: ${e.message ?: e.toString()}
""".stripIndent()
        log.warn("Failed to generate run_stats.txt: ${e.message ?: e.toString()}")
    }
}
