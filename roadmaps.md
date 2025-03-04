# Roadmaps
------
Roadmaps are end-to-end workflows in nextflow. A rule of thumb for a roadmap is to be useful enough on its own to be called seperately. However, it must be designed in a way to work well with other roadmaps. This document provides a list of available roadmaps and their descriptions as well as the instrunctions to run them.

## roadmap_1
### Description
This roadmap is designed to perform metagenomics analysis using a de novo assembly approach. The workflow starts with raw sequencing data and performs the following steps:
1. **Quality Control**: The raw sequencing data is first subjected to quality control using fastp. Next using a reference genome, the reads are decontaminated.
2. **De Novo Assembly**: The cleaned reads are then assembled into contigs using MEGAHIT.
3. **Binning**: The assembled contigs are binned into metagenome-assembled genomes (MAGs) using MetaBAT2.
4. **Abundance Estimation**: The abundance of each bin is estimated using CoverM.

**NOTE**: This roadmap specifically does not include functional or taxonomic annotation of the bins. Those are delegated to other roadmaps.

![roadmap_1](imgs/dag-roadmap_1.svg)

### How to run
Multiple samples can be processed in parallel with nextflow. Roadmap1 workflow needs three inputs:
-   sample_name
-   reads
-   host_genome

Currently there are two ways to run this roadmap:
-   **local**: You have the samples locally in the execution environment. In this case, you need to provide a csv file containing at least three columns:
    -   sample_name
    -   reads1
    -   reads2
Also you need to provide the path to the host genome. An example run with this mode looks like this:
    ```bash
    nextflow run pipelines.nf --roadmap roadmap_1 --host_genome "raw_data/ref_genome.fa" --input_type "local" --input_file <path-to-csv-files> -c configs/local.config
    ```
    **NOTE** You should change the config file according to your environment. 

-   **sra**: In this case you only need a CSV file describig the accession id of your runs:
    -   Run

    An example run with this mode looks like this:
    ```bash
    nextflow run pipelines.nf --roadmap roadmap_1 --host_genome "raw_data/ref_genome.fa" --input_type "sra" --input_file <path-to-csv-files> -c configs/local.config
    ```
