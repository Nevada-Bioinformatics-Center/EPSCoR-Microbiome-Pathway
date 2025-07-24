# EPSCoR NASA Special Project: Pathway-level, consensus analysis of microbiome profiling data

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)


## Introduction

## Pipeline Summary

![](images/pipeline.svg)

## Dataset

### Data 1

Download the data from [**Open Science Data Repository (OSDR)**](https://www.nasa.gov/osdr/). We are currently using the [_OSD-809_](https://osdr.nasa.gov/bio/repo/data/studies/OSD-809) dataset.
Description: Effects of an anaerobic membrance bioreactor upset event on nitrogen speciation and microbial community in a downstream phototrophic membrane bioreactor.

* Factors: 'Time'
* Assay: Metagenomic sequencing - Whole Genome Shotgun Sequencing - Illumina Nextera Kit
* Device Platform: Illumina
* Samples: 12
* Sample Type: Paired-end

### Data 2

Download the data from [**SRA**](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=779554).

_BioProject_: PRJNA779554
Description: Metagenomic analysis of Rhizosphere soil. 
Publication: "A framework for the targeted recruitment of crop-beneficial soil taxa based on network analysis of metagenomic data." doi: 10.1186/s40168-022-01438-1

* Assay: Shotgun Metagenome Sequencing (WGS)
* Device Platform: Illumina NovaSeq 6000
* Samples: 30
* Sample Type: Paired-end

## Installation

> [!NOTE]
> If you are new to Nextflow, please refer to the [this page](https://nextflow.io/docs/latest/install.html) on how to set-up Nextflow.
> We have tested the pipeline on the above mentioned datasets. Please, download the datasets or choose your own dataset.

## Usage

### Basic Usage

1. Run the command below to execute the pipeline.

  ```bash
nextflow main.nf \\
            -profile <PROFILE> \\
            --samplesheet path/to/INPUT.csv \\
            --output /path/to/OUTPUT_DIR \\
            --kneaddata_db /path/to/folder/kneaddata/database/download/ \\
            --metaphlan_db /path/to/folder/metaphlan/database/download \\
            --humann_nucleotide_db /path/to/folder/humann/nucleotide/databases/download \\
            --humann_protein_db /path/to/folder/humann/protein/databases/download \\
            --humann_pathway_db <PATHWAY_DB> \\
            --go_term_db /path/to/folder/go-term/database/download
```

### Parameters

#### Required parameters

| Parameter      | Default              | Description                                                                                              | Options             |
|----------------|----------------------|----------------------------------------------------------------------------------------------------------|---------------------|
| `-profile`     | `local`              | Specify the profile to use for running the pipeline                                                      | `local` , `cluster` |
| `--samplesheet`| NULL                 | Path to a CSV file listing sample names and paired-end FASTQ file paths (R1 and R2), separated by commas |         -           |
| `--output`     | `${baseDir}/results` | Path to the output directory where results will be saved                                                 |         -           |

#### KneadData parameters

| Parameter              | Default | Description                                                                                              | Options |
|------------------------|---------|----------------------------------------------------------------------------------------------------------|---------|
| `--kneaddata_db`       | NULL    | Path to the directory where the KneadData database is located                                            |    -    |
| `--kneaddata_threads`  |    4    | Specify the number of threads                                                                            |    -    |
| `--kneaddata_processes`|    2    | Specify the number of processes                                                                          |    -    |
| `--kneaddata_extra`    |    -    | User specified extra parameter options                                                                   |    -    |

#### Metaphlan parameters

| Parameter           | Default | Description                                                                                              | Options |
|---------------------|---------|----------------------------------------------------------------------------------------------------------|---------|
| `--metaphlan_db`    |    -    | Path to the directory where the MetaPhlAn database is located                                            |    -    |
| `--metaphlan_nproc` |    4    | Specify the number of threads to use                                                                     |    -    |
| `--metaphlan_extra` |    -    | User specified extra parameter options                                                                   |    -    |

#### Humann parameters

| Parameter                | Default   | Description                                                                                              | Options                  |
|--------------------------|-----------|----------------------------------------------------------------------------------------------------------|--------------------------|
| `--humann_nucleotide_db` | NULL      | Path to the directory where the HUMAnN nucleotide database is located                                    |           -              |
| `--humann_protein_db`    | NULL      | Path to the directory where the HUMAnN protein database is located                                       |           -              |
| `--humann_pathway_db`    | `metacyc` | Specify the database to use for pathway computations (default: metacyc)                                  | `metacyc`, `unipathways` |
| `--humann_threads`       |    4      | Specify the number of threads                                                                            |           -              |
| `--humann_extra`         |    -      | User specified extra parameter options                                                                   |           -              |

#### GO-term parameters

| Parameter           | Default | Description                                                                                              | Options |
|---------------------|---------|----------------------------------------------------------------------------------------------------------|---------|
| `--go_term_db`      | NULL    | Path to the directory where the GO-term database is located                                              |    -    |

## Input Files

### Required Files

  1. **Sample Sheet** in `CSV` format. Below are two examples of how the samplesheet should look like:
  
  Example 1: For dataset with exp_conditions (factor) information.

  ```csv
  sample,exp_conditions,fastq_1,fastq_2
  SAMPLE1-ID,NN,sample1_R1.fastq.gz,sample1_R2.fastq.gz
  SAMPLE2-ID,NN,sample2_R1.fastq.gz,sample2_R2.fastq.gz
  SAMPLE3-ID,NN,sample3_R1.fastq.gz,sample3_R2.fastq.gz
  ```

  Example 2: For dataset without exp_conditions (factor) information.

  ```csv
  sample,fastq_1,fastq_2
  SAMPLE1-ID,sample1_R1.fastq.gz,sample1_R2.fastq.gz
  SAMPLE2-ID,sample2_R1.fastq.gz,sample2_R2.fastq.gz
  SAMPLE3-ID,sample3_R1.fastq.gz,sample3_R2.fastq.gz
  ```

  2. **FASTQ files**: Paired-end FASTQ files. In `.fastq.gz` and `.fastq` formats.
   
  3. **KneadData database**: Download the [KneadData database](https://huttenhower.sph.harvard.edu/kneaddata) and provide the path to the database directory.
   
  4. **MetaPhlAn database**: Download the MetaPhlAn database

  5. **HUMAnN nucleotide database**: Download the HUMAnN nucleotide database
   
  6. **HUMAnN protein database**: Download the HUMAnN protein database

## Output


## Citation
