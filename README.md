# PCR Deduplication Analysis Workflow

Evaluates the deduplication process of the `rbt call-consensus-reads` subcommand.
For this evaluation, modifications to to the `call-consensus-reads` source code were made, enabling it to parse name line annotations added to FASTQ files by the ddRAGE software.

## Versions
The release version 1.0.0 is a streamlined version of the code used to generate the plots in Section 7.4.3 of my dissertation.
For easier reproducibility and archival, version 1.1.0 will contain only the parts of rust bio tools actually used in the evaluation, thus reducing unneeded dependencies.

## Requirements

This workflow requires `conda` and `snakemake` to be installed.

For an installation guide for `conda`, please refer to the [conda website](https://docs.conda.io/projects/conda/en/latest/).
This workflow requires the [Python 3 version](https://docs.conda.io/en/latest/miniconda.html) of `(mini)conda`.

We recommend installing `snakemake` via the conda package manager using the bioconda channel to obtain a recent version.
```bash
conda install -c conda-forge -c bioconda snakemake
```
This workflow was tested with `snakemake==5.20.1` on Ubuntu.

The workflow creates approximately 25 GB of intermediate files.

## Usage
Execute the workflow with:

```bash
snakemake --use-conda --jobs <nr cores>
```

The starcode parameters for umi and sequence distances can be configured by editing the file `config.yaml`.
