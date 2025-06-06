nf-metgenomics-piplines: Changelog

The format is based on Keep a Changelog and this project (kind of) adheres to Semantic Versioning.

# Changelog

v1.5.1 - [6/6/25]
* Split to "modules" and "profiles" .configs; defaults to be stored in "modules" and profile-specific values in "profiles"

v1.5.0 - [6/03/25]
* Merged all configs into the modules.config for easier testing and pipeline execution

v1.4.1 - [6/03/25]
* Reformatting to match profile based system

v1.4.0 - [5/30/25]
* Integrating 1 more week of Matt's edits

v1.3.1 - [5/30/25]
* Making tests run with apptainer

v1.3.0 - [5/30/25]
* Integrating branches

v1.2.2 - [5/29/25]
* Allow single_end reads for sylph

v1.2.1 - [5/20/25]
* Added "save_bam" parameter for bowtie2 mapping
* More tests added
* New Docs added for testing and the alternative entry points

v1.2.0 - [5/14/25]
* InStrain subworkflow added
* Centralized configs
* Add "/pl/active/olm-data3/nextflow_scratch" workDir location, and set up automatic deletion after 7 days
* Added "nf-test" for sylph sub_workflow

v1.1.0 - [5/14/25]
* Sylph subworkflow added

v1.0.0 - [5/14/25]
* Added a changelog and versioning to this pipeline
* Added a simple help that displays the version