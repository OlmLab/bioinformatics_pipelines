# nf-metgenomics-piplines

Main repository for  Olm lab bioinformatics pipelines.

Author: Parsa Ghadermazi (pagh1940@colorado.edu)

## Description
-----
This repository contains Nextflow pipelines for computational research in the Olm lab at University of Colorado Boulder. The pipelines are designed to be modular and reusable, allowing for easy integration into various bioinformatics workflows. The main environment to run these pipelines is on Alpine HPC or GutBot internal server.

###  Steps

1. **Install Nextflow**: Make sure you have Nextflow installed on your system. The recommended way to install Nextflow is via conda:
    ```bash
    module load anaconda
    conda create -n nf-pipeline -c bioconda nextflow
    conda activate nf-pipeline
    ```
2. **Clone the repository**: Go to an appropriate location and clone this repository to your local machine using git:
    ```bash
    cd /path/to/your/directory
    git clone https://github.com/ParsaGhadermazi/nf-metgenomics-piplines.git
    ``` 
3. **Pick a roadmap**: Road maps are designed to perform required computations end-to-end. For a full reference, please refer to the [roadmap](roadmaps.md) markdown file.
a roadmap is simply a group of Nextflow workflows chained together. 

4. **Run the pipeline**: To run a pipeline, navigate to the directory containing the pipeline and execute the following command:
    ```bash
    nextflow run pipelines.nf --roadmap <roadmap-id> -c config/alpine.config <roadmap-specific-arguments>
    ```
    Replace `<roadmap-id>` with the name of the pipeline you want to run and `<config_file>` with the path to your configuration file. Each roadmap has its own set of arguments that can be passed to the pipeline. You can find the list of available arguments in the documentation for each [roadmap](roadmaps.md).

or (RECOMENDED) you can skip step 2, and run the pipeline directly by using the following command:

```bash
nextflow run OlmLab/bioinformatics_pipelines --roadmap <roadmap-id> -c config/alpine.config <roadmap-specific-arguments>
```
In this case, the configs will be stored in $HOME/.nextflow/assets/OlmLab/bioinformatics_pipelines/config/

### Optional arguments

- `-profile <profile>`: Specify the execution profile. Most imprtant ones are:
    - `ignore_errors`: Continue execution even if some tasks fail. This is useful for testing and debugging and maybe in some cases when you want skip some of the problematic samples.
- `-resume`: Resume the pipeline from the last completed task. This is useful if the pipeline was interrupted or if you want to re-run only a subset of tasks.


## Recomendations
-----
- Most of the roadmaps require input data in the form of CSV files. Builing such tables can take some time. You can use a tool called bioplumber to make this task easy. Bioplumber is a python package that can be installed via pip:
```bash
pip install bioplumber
```
After installing, you can run the following command to enter the textual user interface for managing your files:
```bash
file-mgr
```
See how this is done as an example 

![bioplumber](/imgs/bioplumber.gif)

- When running pipelines that process many tasks (say 1 million tasks) it is highly recommended that the main nextflow thread is executed from a node with higher memory allocation. Otherwise it might crash!
