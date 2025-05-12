This directory contains the Dockerfiles for the pipelines. Each Dockerfile is designed to create a containerized environment for a specific pipeline or tool used in the Olm lab bioinformatics pipelines. The Dockerfiles are organized by roadmaps. If the dockerfile for a roadmap is not available, it means that the roadmap doesn't need it and relies on other containers for other roadmaps. Here is a brief description of each container:

**NOTES**:

    - The docker files are made in a way to use the latest version of the tools. If you want to use a specific version, you can change the version in the dockerfile.

-------

- **roadmap1**:
    - samtools
    - fastp
    - megahit
    - metabat2
    - SRAToolkit

- **roadmap2**:
    - inStrain
    - samtools
    - bowtie2
    - drep (Only for getting the STB files. THE OTHER FUNCTIONALITIES ARE NOT TESTED)
    - prodigal

- **roadmap3**:
    - drep (Fully Functional)

- **roadmap6**:
    -metaphlan
    -sylph
    -kraken2
    -bracken
    -humann

