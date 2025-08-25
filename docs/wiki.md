# EPSCoR Microbiome Pathway Pipeline Wiki

Welcome to the documentation for the EPSCoR Microbiome Pathway Pipeline. This wiki provides an overview of the pipeline, setup instructions, usage examples, and troubleshooting tips.

## Table of Contents

- [EPSCoR Microbiome Pathway Pipeline Wiki](#epscor-microbiome-pathway-pipeline-wiki)
  - [Table of Contents](#table-of-contents)
  - [Pipeline Overview](#pipeline-overview)
    - [Introduction](#introduction)
    - [Key Features](#key-features)
  - [Getting Started](#getting-started)
    - [System Requirements](#system-requirements)
    - [Installing Nextflow](#installing-nextflow)
    - [Tools and Dependencies](#tools-and-dependencies)
  - [Usage](#usage)
  - [Pipeline Steps](#pipeline-steps)
    - [Kneading Data](#kneading-data)
    - [Taxonomic Profiling](#taxonomic-profiling)
    - [Merge Taxonomic Tables](#merge-taxonomic-tables)
    - [Taxonomic Visualization](#taxonomic-visualization)
    - [Functional Profiling](#functional-profiling)
    - [Normalize Pathway Abundance](#normalize-pathway-abundance)
    - [Join Pathway Abundance](#join-pathway-abundance)
    - [Descriptive Profiling](#descriptive-profiling)
    - [Generate GOterms](#generate-goterms)
    - [Consensus Pathway Analysis](#consensus-pathway-analysis)
    - [Run MultiQC](#run-multiqc)
  - [Configuration](#configuration)
  - [Helpful Tips](#helpful-tips)

## Pipeline Overview

### Introduction

*{Pipeline-name}* is a scalable and reproducible pipeline for comprehensive taxonomic, functional, and consensus pathway analyses of metagenomic data. Built using [Nextflow](https://nextflow.io), the pipeline supports execution on local machines, as well as high-performance cluster (HPC) and cloud environments. Several popular bioinformatics tools and packages are integrated in the pipeline to perform tasks such as sequence quality control, taxonomic, and functional profiling, and downstream statistical analyses.

### Key Features

- Performs consensus pathway enrichment analyses on metagenomic data
- Descriptive profiling of pathways, respective taxa, and abundances present in samples
- Alpha and Beta diversity analysis by generating publication-ready plots
- Taxonomic and Functional Profiling
- Comprehensive QC reports

## Getting Started

### System Requirements

The pipeline is currently capable of running on any UNIX-based system, such as Linux, macOS, etc. For running the pipeline locally, it is recommended to use a computer with at least 64GB of memory and a storage of *>=1 TB*. To run the pipeline on HPC or cloud systems such as AWS, there is no recommended specification.

Since the pipeline is built using Nextflow, the user must install it. Nextflow installation can be found in the next sub-section. Additionally, the user must also ensure to install [Conda](https://docs.conda.io/en/latest/) when running locally and [Docker](https://www.docker.com) when running on other platforms.

### Installing Nextflow

Before installing Nextflow, please ensure that **Bash** version *3.2 (or later)* and **Java** version *17 (or later, up to 24)* are installed on your system. If these are not present, you will need to install them.

It is recommended to install Nextflow using the self-installing package distribution mode; however, users can also opt for installation via conda mode or as a standalone distribution. The following steps outline the process for installing Nextflow as a self-installing package:

**Step 1:** Download Nextflow

```bash
curl -s https://get.nextflow.io | bash
```

**Step 2:** Make Nextflow executable

```bash
chmod +x nextflow
```

**Step 3:** Move Nextflow to an executable path. For example -

```bash
mkdir -p $HOME/.local/bin/
mv nextflow $HOME/.local/bin/
```

**Step 4:** Finally, verify that Nextflow is installed correctly

```bash
nextflow info
```

The pipeline utilizes **25.04.6 build 5954**. For more information, please visit: https://nextflow.io/docs/latest/install.html

### Tools and Dependencies

The pipeline uses widely recognized bioinformatics tools and packages to conduct various tasks, including sequence quality control, taxonomic and functional profiling, as well as downstream statistical analyses and visualizations.

Currently, the pipeline incorporates the following tools and R packages:

- KneadData (version 0.12.2)
- FastQC (version 0.12.1)
- MetaPhlAn (version 4.1.1)
- HUMAnN (version 3.9)
- MultiQC (version 1.29)
- R Base
- Pacman
- dplyr
- ggplot2
- Bioconductor/BiocManager
- vegan
- mia
- phyloseq
- DirichletMultinomial
- stringr
- Tidyverse
- fgsea
- limma
- SummarizedExperiment
- Uniprot.ws
- Devtools
- igraph
- ggraph
- RColorBrewer

> [!NOTE]
> All tools and dependencies are installed automatically by the pipleine using Conda and YAML files when running locally, and via Docker containers when executed on other platforms.
> These files can be found in the `./assets` directory.

## Usage

Before executing the pipeline, please follow these mandatory steps:

**Step 1:** Prepare a sample sheet in CSV format that lists the samples and their metadata. Here are two examples of the sample sheet:

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

**Step 2:** Paired-end FASTQ files in either `.fastq` or `.fastq.gz` format.

**Step 3:** Provide a database for host sequence removal. This can be downloaded from the KneadData database by running the following command:

```bash
kneaddata_database --download <DATABASE> <BUILD> <DATABASE_FOLDER>
```

To view the list of available databases, use:

```bash
kneaddata_database --available
```

Alternatively, the user can build/provide custom reference database. For more information, see the [KneadData documentation](https://github.com/biobakery/kneaddata)

**Step 4:** Download the MetaPhlAn database by running the following command:

```bash
metaphlan --install --index <INDEX> --bowtie2db <DATABASE_FOLDER>
```

The default index is set to `mpa_vJun23_CHOCOPhlAnSGB_202403`. Alternatively, the user may download the database from the [Segata Lab FTP site](https://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/?C=M;O=D). When donwloading from the FTP site, please ensure to untar the files.

> [!WARNING]
> Do not use the latest MetaPhlAn database, as it is incompatible with HUMAnN 3.9.
> For further information, consult the [BioBakery forum](https://forum.biobakery.org).

**Step 5:** Download both the HUMAnN nucleotide and protein databases using the following command:

```bash
humann_databases --download <DATABASE> <BUILD> <DIRECTORY>
```

To view all available databases, use:

```bash
humann_databases --available
```

For more information, please review the [HUMAnN documentation](https://github.com/biobakery/humann?tab=readme-ov-file#5-download-the-databases).

**Step 6:** Finally, execute the pipeline using the following command:

```bash
nextflow main.nf \\
-profile <PROFILE> \\
--samplesheet <PATH-to-INPUT.csv> \\
--output <PATH-to-OUTPUT_DIR> \\
--kneaddata_db <PATH-to-folder-KneadData_database> \\
--metaphlan_db <PATH-to-folder-MetaPhlAn_database> \\
--humann_nucleotide_db <PATH-to-folder-HUMAnN_nucleotide_database> \\
--humann_protein_db <PATH-to-folder-HUMAnN_protein_database> \\
--humann_pathway_db <PATHWAY_DB> \\
--go_term_db <PATH-to-folder-GOTerm_database
```

## Pipeline Steps

### Kneading Data

The Kneading Data process performs comprehensive quality control and host decontamination of raw metagenomic sequence reads. This step ensures that only high-quality, host-depleted reads are retained for downstream analyses, thus enhancing the accuracy of taxonomic and functional profiling.

*Key Features:*

- Trims low-quality bases and removes adapters using Trimmomatic via KneadData.
- Removes host and contaminant sequences by aligning reads to reference genomes with Bowtie2.
- Accepts paired-end FASTQ files and outputs merged, host-depleted FASTQ files.
- Runs FastQC before and after processing for quality assessment.
- Supports custom reference databases, sequencer source, thread count, and KneadData options.
- Cleans up intermediate files and organizes outputs for downstream analysis.

*Parameters:*

| Parameter                        | Description                                                                 | Default    |
|-----------------------------------|-----------------------------------------------------------------------------|------------|
| `--kneaddata_db`                  | Path to the KneadData database directory (required).                        | None       |
| `--kneaddata_sequencer_source`    | Sequencer type: NexteraPE, TruSeq2, TruSeq3 or none                         | NexteraPE  |
| `--kneaddata_threads`             | Number of threads to use per process.                                       | 4          |
| `--kneaddata_processes`           | Number of parallel processes to run.                                        | 1          |
| `--kneaddata_time`                | Maximum time allowed for the process (used only when cluster profile is selected) | 8h         |
| `--kneaddata_mem`                 | Maximum memory allowed for the process.                                     | 8 GB       |
| `--kneaddata_extra`               | Additional user-specified options to pass directly to KneadData.            | (Empty)    |

*Outputs:*

- `${sample}.fastq.gz` : Host-depleted, quality-controlled FASTQ file for each sample.
- `${sample_id}_kneaddata.log` : Detailed logs for each run.
- `.html` and `.zip` reports for processed reads, generated by FastQC after KneadData processing.
- `kneaddata_out/` : Directory containing all KneadData output files organized by sample.

> [!NOTE]
> All the outputs have symbolic links in the `./results/kneaddata/` directory. 
> If you wish to access the actual files, please navigate to the `work/` directory and locate the respective folders.

### Taxonomic Profiling

The Taxonomic Profiling process identifies and quantifies the microbial taxa present in metagenomic samples using MetaPhlAn. This step provides detailed insights into the composition of microbial communities at various taxonomic levels, such as kingdom, phylum, class, order, family, genus, and species.

*Key Features:*

- Accurate taxonomic assignment using MetaPhlAn’s marker gene approach for high-resolution, low-false-positive profiles.
- Outputs relative abundance tables, Bowtie2 alignment files, and concise log files for each sample.
- Supports multi-threading, process parallelization, and customizable MetaPhlAn options.

*Parameters:*

| Parameter                  | Description                                                                 | Default    |
|----------------------------|-----------------------------------------------------------------------------|------------|
| `--metaphlan_db`           | Path to the MetaPhlAn database directory (required).                        | None       |
| `--metaphlan_index`        | MetaPhlAn database index to use.                                           | mpa_vJun23_CHOCOPhlAnSGB_202403 |
| `--metaphlan_threads`      | Number of threads to use per process.                                       | 4          |
| `--metaphlan_processes`    | Number of parallel processes to run.                                        | 1          |
| `--metaphlan_time`         | Maximum time allowed for the process (used only when cluster profile is selected) | 8h         |
| `--metaphlan_mem`          | Maximum memory allowed for the process.                                     | 16 GB       |
| `--metaphlan_extra`        | Additional user-specified options to pass directly to MetaPhlAn.            | (Empty)    |

*Outputs:*

- `${sample_id}_profile.tsv` and `.txt` : Taxonomic table containing relative abundances of identified taxa.
- `${sample_id}_metaphlan.log` : Log file summarizing only warnings and errors, else empty.
- `metaphlan_out/` : Directory containing all MetaPhlAn output files organized by sample.

### Merge Taxonomic Tables

The Merge Taxonomic Tables process combines individual taxonomic abundance tables generated by MetaPhlAn for each sample into a single, comprehensive table. This merged table provides a unified view of taxonomic profiles across all samples, facilitating downstream comparative analyses and visualization.

*Key Features:*

- Uses MetaPhlAn’s `merge_metaphlan_tables.py` to efficiently combine taxonomic profiles from all samples.
- Includes all detected taxa, filling missing values as needed.
- Outputs a standardized, analysis-ready table for downstream visualization and statistics.

*Outputs:*

- `merged_taxa_profile.tsv` : A single file containing the relative abundances of all detected taxa across all samples. Each row represents a taxon, and each column corresponds to a sample.
- `metaphlan_out/merge` : Directory containing the merged taxonomic profile.

### Taxonomic Visualization

The Taxonomic Visualization process generates a comprehensive set of plots and summary statistics from the merged taxonomic abundance table produced by MetaPhlAn. This step enables users to visually explore and interpret the microbial community composition across all samples.

*Key Features:*

- Generates publication-ready plots (e.g., stacked bar plots at phylum/genus levels, alpha diversity, NMDS ordination) using an R script.
- Accepts the merged taxonomic profile and samplesheet for grouping and comparison by experimental conditions.
- Produces figures for interpretation and reporting, supporting user-defined grouping and visualization options.
- Saves all visualizations in a dedicated output directory for easy access.

*Outputs:*

- `metaphlan_out/plots/` : Directory containing all generated plots in `.png` format.

> [!NOTE]
> To enable both the Merge Taxonomic Tables and Taxonomic Visualization processes, set the boolean parameter `--metaphlan_extra_analysis` to `true` when executing the pipeline.

### Functional Profiling

The Functional Profiling process characterizes the functional potential of metagenomic samples using HUMAnN. This step identifies and quantifies gene families and metabolic pathways present in the microbial community, providing insights into the biological functions encoded by the microbiome.

*Key Features:*

- Functional profiling with HUMAnN, mapping reads to gene families and reconstructing pathways using nucleotide and protein databases.
- Integrates MetaPhlAn taxonomic profiles for improved accuracy.
- Outputs gene family, pathway abundance, and coverage tables per sample.
- Supports user-specified databases (nucleotide, protein, pathway).
- Enables multi-threading and process parallelization.
- Flexible parameters for threads, database selection, and extra HUMAnN options.
- Organizes all results and logs in a dedicated output directory.

*Parameters:*

| Parameter                  | Description                                                                 | Default    |
|----------------------------|-----------------------------------------------------------------------------|------------|
| `--humann_nucleotide_db`   | Path to the HUMAnN nucleotide database directory (required).                | None       |
| `--humann_protein_db`      | Path to the HUMAnN protein database directory (required).                   | None       |
| `--humann_pathway_db`      | Pathway database to use: `metacyc` or `unipathways`.                       | metacyc    |
| `--humann_threads`         | Number of threads to use per process.                                       | 4          |
| `--humann_time`            | Maximum time allowed for the process (used only when cluster profile is selected). | 8h         |
| `--humann_mem`             | Maximum memory allowed for the process.                                     | 16 GB      |
| `--humann_extra`           | Additional user-specified options to pass directly to HUMAnN.               | (Empty)    |

*Outputs:*

- `${sample_id}_genefamilies.tsv` : Gene family abundance table for each sample.
- `${sample_id}_pathabundance.tsv` : Pathway abundance table for each sample.
- `${sample_id}_pathcoverage.tsv` : Pathway coverage table for each sample.
- `${sample_id}_humann.log` : Log file summarizing details of the HUMAnN run.
- `humann_out/` : Directory containing all HUMAnN output files organized by sample.

### Normalize Pathway Abundance

The Normalize Pathway Abundance process standardizes the pathway abundance tables generated by HUMAnN for each sample. This normalization step converts raw pathway counts into relative abundances or counts per million (CPM), enabling meaningful comparisons across samples.

*Key Features:*

- Utilizes `humann_renorm_table` to convert pathway abundances from one normalization method to another.
- Supports output in relative abundance (relab) or counts per million (cpm).
- Saves normalized tables in a dedicated output directory per sample.

*Parameters:*

| Parameter              | Description                                      | Default |
|------------------------|--------------------------------------------------|---------|
| `--humann_renorm_units`| Units for normalization: `relab` or `cpm`        | cpm     |

*Outputs:*

- `${sample_id}_renorm_pathabundance.tsv` : Normalized pathway abundance table for each sample.
- `humann_out/pathabundance` : Directory containing all normalized pathway abundance tables organized by sample.

### Join Pathway Abundance

The Join Pathway Abundance process merges the normalized pathway abundance tables from all samples into a single, comprehensive matrix. This joined table facilitates cohort-level analyses and downstream statistical comparisons.

*Key Features:*

- Automated merging of normalized pathway tables using `humann_join_tables`.
- Produces a unified matrix with all pathways across samples, handling missing values.
- Output is ready for downstream analysis and visualization.

*Outputs:*

- `joined_norm_pathabundance.tsv` : A single file containing the normalized abundances of all detected pathways across all samples. Each row represents a pathway, and each column corresponds to a sample.
- `humann_out/merge` : Directory containing the joined normalized pathway abundance table.

### Descriptive Profiling

The Descriptive Profiling process summarizes the joined pathway abundance data, providing insights into the most abundant taxa and pathways across the dataset. This step helps users quickly identify key functional features and their taxonomic contributors.

*Key Features:*

- Automated R-based summarization of joined pathway abundance data.
- Identifies top genus-pathway associations by maximum abundance.
- Outputs organized results for easy review and downstream analysis.

*Outputs:*

- `top_path_taxa_results.tsv` : A summary table listing the top genus-pathway associations based on maximum abundance across all samples.
- `desc_out/` : Directory containing the descriptive profiling result file.

### Generate GOterms

The Generate GO Terms process prepares the gene ontology (GO) resources required for consensus pathway analysis (CPA) of metagenomic data. This step downloads, parses, and formats the necessary annotation files, producing a ready-to-use gene set object for downstream functional enrichment analysis.

*Key Features:*

- Automated download of gene2go, gene_info, and GO OBO files from NCBI and the Gene Ontology Consortium.
- Efficient parsing to extract gene-to-GO mappings and biological process terms.
- Construction of gene sets for enrichment analysis, filtering by size.
- Outputs all files and a serialized R object (`GOTerms.rds`) for downstream use.
- Skips steps if files already exist for efficient re-runs.

*Parameters:*

| Parameter                   | Description                                                                 | Default    |
|-----------------------------|-----------------------------------------------------------------------------|------------|
| `--goterm_db`               | Path to the directory where GO term files will be downloaded and stored.     | None       |
| `--goterm_cpu`              | Number of CPUs to use for GO term processing.                               | 2          |
| `--goterm_mem`              | Maximum memory allowed for GO term processing.                              | 16 GB      |

*Outputs:*

- `GOTerms.rds` : Serialized R object containing the gene sets and GO term names for CPA.
- `gene2go.gz` : Downloaded gene-to-GO mapping file from NCBI.
- `All_Data.gene_info.gz` : Downloaded gene information file from NCBI.
- `go.obo` : Downloaded ontology structure file from the Gene Ontology Consortium.
- `.done` : Marker file indicating successful completion of the GO term generation process.

### Consensus Pathway Analysis

The Consensus Pathway Analysis (CPA) process integrates gene family abundance data with gene ontology (GO) resources to identify and visualize functionally significant pathways across experimental conditions in metagenomic datasets. This step performs differential expression and gene set enrichment analyses, providing robust statistical and visual summaries of pathway-level changes.

*Key Features:*

- Differential expression analysis of gene family abundances between groups (e.g., using limma)
- Gene set enrichment analysis for GO terms via multiple methods (ORA, FGSEA, KS, Wilcoxon)
- Consensus scoring to highlight the most significant pathways
- Publication-ready network and bubble plots for top GO terms and pathways
- All results and figures organized in a dedicated output directory

*Parameters:*

| Parameter     | Description                                                                                   | Default |
|---------------|-----------------------------------------------------------------------------------------------|---------|
| `--cpa_cpu`   | Number of CPUs to use for Consensus Pathway Analysis processing.                              | 2       |
| `--cpa_mem`   | Maximum memory allowed for Consensus Pathway Analysis processing.                             | 16 GB   |
| `--cpa_time`  | Maximum time allowed for the Consensus Pathway Analysis process (used with cluster profiles). | 24h     |

*Outputs:*

- `GOTerms.rds` : Serialized R object containing the gene sets and GO term names for consensus pathway analysis.
- `cpa_results.csv` : CSV file summarizing pathway enrichment and differential expression results across experimental conditions.
- `network_plot.png` : Publication-ready network plot visualizing the top GO terms and their relationships.
- `bubble_plot.png` : Bubble plot highlighting the most significant GO terms across comparisons.
- `cpa_out/` : Directory containing all CPA results, figures, and logs organized for easy access.

### Run MultiQC

The MultiQC process aggregates quality control (QC) metrics and summary statistics from all upstream steps in the pipeline, including raw and processed read QC reports, into a single, interactive report. This step provides a comprehensive overview of data quality and processing outcomes across all samples in the analysis.

*Key Features:*

- Aggregates QC reports (FastQC, KneadData, Trimmomatic) into a single interactive HTML summary.
- Automatically detects and includes all relevant QC outputs from raw and processed data.
- Saves the MultiQC report and data in a dedicated output directory for easy access.

*Outputs:*

- `multiqc_report.html`: An interactive HTML report summarizing QC metrics for all samples.
- `multiqc_data/`: Directory containing the raw data files used to generate the MultiQC report.
- `multiqc_out/`: Directory containing the MultiQC report and associated data files.

## Configuration

## Helpful Tips

- Nextflow has the ability to cache task executions and re-use them in later runs so to minimize duplicate work. If encountered error or for any other reason, the user can resume the pipeline with the `-resume` flag.

  ```bash
  nextflow run main.nf -resume
  ```
- The pipeline will not automatically delete the work directory. The user may choose to either keep the directory and clean the cache by the following command:

```bash
nextflow clean [run_name|session_id] [options]
```
or, permanently delete it by:

```bash
rm -rf work
```
For more information, please see [here](https://www.nextflow.io/docs/latest/reference/cli.html#clean).