# nf-metgenomics-piplines

Main repository for the Olm lab bioinformatics pipelines.

Author: Parsa Ghadermazi (pagh1940@colorado.edu)

## Description
-----
This repository contains Nextflow pipelines for computational research in the Olm lab at University of Colorado Boulder. The pipelines are designed to be modular and reusable, allowing for easy integration into various bioinformatics workflows. The main environment to run these pipelines is on Alpine HPC or GutBot internal server.

###  Steps

1. **Install Nextflow**: Make sure you have Nextflow installed on your system. The recommended way to install Nextflow is via conda:
    ```bash
    moduel load anaconda
    conda create -n nf-pipeline -c bioconda nextflow
    conda activate nf-pipeline
    ```
2. **Clone the repository**: Go to an appropriate location and clone this repository to your local machine using git:
    ```bash
    cd /path/to/your/directory
    https://github.com/OlmLab/bioinformatics_pipelines.git
    ``` 
3. **Pick a roadmap**: Road maps are designed to perform required computations end-to-end. For a full reference, please refer to the [roadmap]() markdown file.
a roadmap is simply a group of Nextflow workflows chained together. 

4. **Run the pipeline**: To run a pipeline, navigate to the directory containing the pipeline and execute the following command:
    ```bash
    nextflow run pipelines.nf --roadmap <roadmap-id> -c config/alpine.config <roadmap-specific-arguments>
    ```
    Replace `<roadmap-id>` with the name of the pipeline you want to run and `<config_file>` with the path to your configuration file. Each roadmap has its own set of arguments that can be passed to the pipeline. You can find the list of available arguments in the documentation for each [roadmap]().
