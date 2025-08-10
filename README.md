# EPSCoR NASA Special Project: Pathway-level, consensus analysis of microbiome profiling data

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)


## Introduction

## Pipeline Summary

![](images/pipeline.png)

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

#### Global Options

| Parameter         | Default                | Description                                                                                  | Options                |
|-------------------|------------------------|----------------------------------------------------------------------------------------------|------------------------|
| `-profile`        | `local`                | Execution profile                                                                            | `local`, `cluster`     |
| `--samplesheet`   | NULL                   | CSV file with sample names and paired-end FASTQ file paths (R1, R2)                          | -                      |
| `--output`        | `${baseDir}/results`   | Output directory for results                                                                 | -                      |

#### KneadData Parameters

| Parameter                   | Default      | Description                                                      | Options                        |
|-----------------------------|--------------|------------------------------------------------------------------|--------------------------------|
| `--kneaddata_db`            | NULL         | Path to KneadData database                                       | -                              |
| `--kneaddata_sequencer_source` | NexteraPE | Sequencer type                                                   | `NexteraPE`, `TruSeq2`, `TruSeq3` |
| `--kneaddata_threads`       | 4            | Threads for KneadData                                            | -                              |
| `--kneaddata_processes`     | 1            | Processes for KneadData                                          | -                              |
| `--kneaddata_time`          | -            | Process time (cluster only)                                      | -                              |
| `--kneaddata_mem`           | 8 GB         | Memory for KneadData                                             | -                              |
| `--kneaddata_extra`         | -            | Additional KneadData options                                     | -                              |

#### Taxonomic Profiling (MetaPhlAn) Parameters

| Parameter                   | Default                              | Description                                                      | Options                        |
|-----------------------------|--------------------------------------|------------------------------------------------------------------|--------------------------------|
| `--metaphlan_db`            | -                                    | Directory for MetaPhlAn database                                 | -                              |
| `--metaphlan_index`         | mpa_vJun23_CHOCOPhlAnSGB_202403      | MetaPhlAn index                                                  | -                              |
| `--metaphlan_nproc`         | 4                                    | Processes for MetaPhlAn                                          | -                              |
| `--metaphlan_time`          | -                                    | Process time (cluster only)                                      | -                              |
| `--metaphlan_mem`           | 16 GB                                | Memory for MetaPhlAn                                             | -                              |
| `--metaphlan_extra`         | -                                    | Additional MetaPhlAn options                                     | -                              |
| `--metaphlan_extra_analysis`| -                                    | Estimate Shannon, Inverse Simpson, and Alpha diversity           | -                              |

#### Functional Profiling (HUMAnN) Parameters

| Parameter                   | Default      | Description                                                      | Options                        |
|-----------------------------|--------------|------------------------------------------------------------------|--------------------------------|
| `--humann_nucleotide_db`    | NULL         | HUMAnN nucleotide database directory                             | -                              |
| `--humann_protein_db`       | NULL         | HUMAnN protein database directory                                | -                              |
| `--humann_pathway_db`       | metacyc      | Pathway DB for HUMAnN                                            | `metacyc`, `unipathways`       |
| `--humann_threads`          | 4            | Threads for HUMAnN                                               | -                              |
| `--humann_time`             | -            | Process time (cluster only)                                      | -                              |
| `--humann_mem`              | 16 GB        | Memory for HUMAnN                                                | -                              |
| `--humann_extra`            | -            | Additional HUMAnN options                                        | -                              |

#### Consensus Pathway Analysis Parameters

| Parameter                   | Default      | Description                                                      | Options                        |
|-----------------------------|--------------|------------------------------------------------------------------|--------------------------------|
| `--goterm_db`               | NULL         | GOterms database directory                                       | -                              |
| `--goterm_cpu`              | 2            | CPUs for GOterm analysis                                         | -                              |
| `--goterm_mem`              | 16 GB        | Memory for GOterm analysis                                       | -                              |
| `--cpa_cpu`                 | 2            | CPUs for consensus pathway analysis                              | -                              |
| `--cpa_mem`                 | 16 GB        | Memory for consensus pathway analysis                            | -                              |

#### Descriptive Profiling Analysis Parameters

| Parameter                   | Default      | Description                                                      | Options                        |
|-----------------------------|--------------|------------------------------------------------------------------|--------------------------------|
| `--humann_renorm_units`     | cpm          | Pathway abundance normalization units                            | `cpm`, `relab`                 |
| `--norm_path_cpu`           | 2            | CPUs for descriptive profiling                                   | -                              |
| `--norm_path_mem`           | 16 GB        | Memory for descriptive profiling                                 | -                              |
| `--join_path_cpu`           | 1            | CPUs for join pathway                                            | -                              |
| `--join_path_mem`           | 8 GB         | Memory for join pathway                                          | -                              |
| `--desc_cpu`                | 1            | CPUs for descriptive analysis                                    | -                              |
| `--desc_mem`                | 8 GB         | Memory for descriptive analysis                                  | -                              |

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
