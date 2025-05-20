# Entries
------
Entry points are an alternative way of accessing the modules within this codebase. You can get a view of all the current entry points by running the help:

```
$ nextflow main.nf --help

 N E X T F L O W   ~  version 24.10.6

Launching `main.nf` [mighty_faraday] DSL2 - revision: 646250ea80


Olm Lab nf-metgenomics-pipelines Pipeline Version: 1.2.1

Usage: nextflow run main.nf --entry <subworkflow or workflow> [options]

Available Workflows:

  - wf_instrain_profile

Available Subworkflows:

  - sw_sylph
  - sw_bowtie2
```

You can then run a particular workflow using a command like:

```
nextflow run /pl/active/olm-data1/software/olm/bioinformatics_pipelines/main.nf -entry wf_instrain_profile --input_file NeMaSiMa_passScreen_v1.csv  -resume --genome_db /pl/active/olm-data2/Data/Collab/TinyHealth2025/phase1/other/output_bifido_gtdb/genomes_db.fasta --stb /pl/active/olm-data2/Data/Collab/TinyHealth2025/phase1/other/output_bifido_gtdb/instrain/stb/genomes_db.stb -c /pl/active/olm-data1/software/olm/bioinformatics_pipelines/configs/blanca.config --output_dir NeMaSiMa_profiles_v1
```

The availble input parameters for various workflows and subworkflows can be found in the main [nextflow config file](../nextflow.config)