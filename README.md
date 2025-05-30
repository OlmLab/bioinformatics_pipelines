# nf-metgenomics-piplines

Main repository for  Olm lab bioinformatics pipelines.

Author: Parsa Ghadermazi (pagh1940@colorado.edu)

## Description
-----
This repository contains Nextflow pipelines for computational research in the Olm lab at University of Colorado Boulder. The pipelines are designed to be modular and reusable, allowing for easy integration into various bioinformatics workflows. 

There are currently two ways to run this code- one of which is based on roadmaps, and one of which is based on entry points. They have different documentation, as outlined below.

[Documentation for roadmap-based workflows](docs/roadmaps.md)

[Documentation for entry-based workflows](docs/entries.md)

## Testing
-----

Basic testing is provided using the [nf-test](https://www.nf-test.com) framework. Run the following command from the base directory of this repo to run all tests:

```
$ nf-test test
```


## Bioplumber
-----

Most of the roadmaps require input data in the form of CSV files. Builing such tables can take some time. You can use a tool called bioplumber to make this task easy. Bioplumber is a python package that can be installed via pip:

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
