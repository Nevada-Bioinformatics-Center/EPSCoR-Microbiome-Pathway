# EPSCoR NASA Special Project: Pathway-level, consensus analysis of microbiome profiling data

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)


## Introduction

## Pipeline Summary

![](images/pipeline.png)

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

**Step 1:** Prepare a **Sample Sheet (CSV format)** listing your samples and metadata information.

 With experimental conditions (factors):

```csv
      sample,exp_conditions,fastq_1,fastq_2
      SAMPLE1-ID,NN,sample1_R1.fastq.gz,sample1_R2.fastq.gz
      SAMPLE2-ID,NN,sample2_R1.fastq.gz,sample2_R2.fastq.gz
      SAMPLE3-ID,NN,sample3_R1.fastq.gz,sample3_R2.fastq.gz
```

 Without experimental conditions:

```csv
      sample,fastq_1,fastq_2
      SAMPLE1-ID,sample1_R1.fastq.gz,sample1_R2.fastq.gz
      SAMPLE2-ID,sample2_R1.fastq.gz,sample2_R2.fastq.gz
      SAMPLE3-ID,sample3_R1.fastq.gz,sample3_R2.fastq.gz
```

**Step 2:** Provide **Paired-end FASTQ files** in `.fastq.gz` or `.fastq` format.

**Step 3:** Download **KneadData Database** by running the following command

```bash
    kneaddata_database --download <DATABASE> <BUILD> <DATABASE_FOLDER>
```

To view the list of available databases:

```bash
    kneaddata_database --available
```

Alternatively, you can also build your own custom reference database. See [KneadData README](https://github.com/biobakery/kneaddata) for more information.

> [!NOTE]
> Only Bowtie2-generated databases are supported.

**Step 4:** Download **MetaPhlAn Database** by the following command

```bash
    metaphlan --install --index <INDEX> --bowtie2db <DATABASE_FOLDER>
```

Presently, the default index has been set to `mpa_vJun23_CHOCOPhlAnSGB_202403`. 

Alternatively, you can also download from [Segata Lab FTP](http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/?C=M;O=D). If following this, please do untar the file.

> [!CAUTION]
> Do **not** use the latest MetaPhlAn database; it is incompatible with HUMAnN 3.9. 
> For more information please consult the [BioBakery forum](https://forum.biobakery.org/)

For more information please see [MetaPhlAn README](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-4.1).

**Step 5:** Download both **HUMAnN Nucleotide and Protein Databases** by the following command

```bash
    humann_databases --download <DATABASE> <BUILD> <DIRECTORY>
```

To view the available databases:

```bash
    humann_databases --available
```

For more information, please review [HUMAnN README](https://github.com/biobakery/humann?tab=readme-ov-file#5-download-the-databases)

## Output

The pipeline generates the following output files:

> [!NOTE]
> The output files are stored in nested directories within the user-defined output directory.

1. **kneaddata_out**
    -  Processed (decontaminated, trimmed) FASTQ files in `.fastq.gz` format. 
    -  `_kneaddata.log` file containing the log of the KNEADING_DATA process for each samples.
    -  `/fastqc` directory containing FastQC reports for the processed FASTQ files in `.html` and the corresponding `.zip` files.

> Please note that these have symbolic links to the nextflow working directory.

2. **metaphlan_out**
    -  `_profile.tsv` file containing the taxonomic profile of the samples.
    -  `_profile.txt` file containing the taxonomic profile in text format.
    -  `_metaphlan.log` file will contain warnings and errors from the MetaPhlAn process. Else it is empty.
    -  `/merge` directory containing the merged taxonomic profile across all samples in `_profile.tsv` format. This file is in `tsv` format.
    -  `/plots` directory containing the taxonomic profile plots in `.png` format. The plots include:
        -  `genus_abundance.png`: Genus-level abundance plot.
        -  `phyla_abundance.png`: Phyla-level abundance plot.
        -  `scaled_abundances_top10.png`: Scaled abundance plot for the top 10 taxa.
        -  `Bray_NMDS.png`: Bray-Curtis NMDS plot.
        -  `alpha_diversity_shannon.png`: Shannon diversity plot.
        -  `NMDS.png`: NMDS plot of the taxonomic profile.
        -  `alpha_diversity_inverse_simpson.png`:  Inverse Simpson diversity plot

3. **humann_out**
    -  `_genefamilies.tsv` file containing the gene families profile.
    -  `_pathabundance.tsv` file containing the pathway abundance profile.
    -  `_pathcoverage.tsv` file containing the pathway coverage profile.
    -  `_humann.log` file contains the log of the FUNCTIONAL_PROFILING process for each sample.
    -  `/pathabundance` directory containing the counts per million (CPM) converted pathway abundance profile across all samples in `_renorm_pathabundance.tsv` format.
    -  `/merge` directory containing the merged CPM converted pathway abundance profile across all samples in `.tsv` format.

4. **cpa_out**
    -  `abundances.csv` file containing the gene families abundances across all samples.
    -  `mapping.csv` file containing the mapping of gene families to UniProt.
    -  `mappingAbundanceMatrix.csv` file containing the mapping of gene families to UniProt with abundance values.
    -  `goterms-results.csv` file containing the GO terms results with ID, p-value, normalized score, pFDR, and description.
    -  `consensus-results.csv` file containing the consensus results with ID, description of pathways, and p-fisher value.
    -  `/export_Knead` directory containing the pair-wise comparison of experimental conditions. The files are in `.csv` format.
    -  `Pathway_Bubble_Plot.png` file containing the bubble plot of the pathways.
    -  `Pathway_Networks_Plot.png` file containing the network plot of the pathways.

5. **desc_out**
    - `top_path_taxa_results.tsv` file containing the top pathway taxa results.

6. **multiqc_out**
   - `multiqc_report.html` file containing the MultiQC report of the pipeline run.
   - `/multiqc_data` directory containing the MultiQC data files.

## Citation
