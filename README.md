# EPSCoR NASA Special Project: Pathway-level, consensus analysis of microbiome profiling data

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-2496ED?labelColor=000000&logo=docker)](https://docs.docker.com/get-docker/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1F1F1F?labelColor=000000&logo=singularity)](https://sylabs.io/guides/)


## Introduction

This is a scalable and reproducible pipeline, built in [Nextflow](https://www.nextflow.io/) for comprehensive taxonomic, functional, and consensus pathway analysis of metagenomes. It supports execution on local machines, as well as HPC and cloud environments. It integrates established tools for sequence quality control, taxonomic profiling, and functional profiling, as well as a consensus-based pathway enrichment strategy that combines multiple statistical approaches to improve robustness. In addition, the pipeline compiles a descriptive summary that links microbial taxa with associated functions.

## Pipeline Summary

![](images/pipeline.png)

## Usage

> [!TIP]
> For detailed information on installation, usage, outputs & troubleshooting, please refer to the [wiki document](docs/wiki.md).

> [!NOTE]
> If you are new to Nextflow, please refer to the [this page](https://nextflow.io/docs/latest/install.html) on how to set-up Nextflow.

### Basic Usage

Run the command below to execute the pipeline:

```bash
nextflow run main.nf \
-profile <EXECUTION-PROFILE>,<ENGINE-TOGGLE> \
--samplesheet path/to/INPUT.csv \
--output /path/to/OUTPUT_DIR \
--kneaddata_db /path/to/folder/kneaddata/database/download/ \
--metaphlan_db /path/to/folder/metaphlan/database/download \
--humann_nucleotide_db /path/to/folder/humann/nucleotide/databases/download \
--humann_protein_db /path/to/folder/humann/protein/databases/download \
--humann_pathway_db <PATHWAY_DB> \
--goterm_db /path/to/folder/go-term/database/download \
-resume
```

For help on the available pipeline parameters, run the command below:

```bash
nextflow main.nf --help
```

For the impatient, real usage:

```bash
nextflow run main.nf \
-profile slurm,singularity \
--samplesheet SAMPLESHEET.csv \
--output results \
--kneaddata_db /home/mypool/projects/nasa_pipeline/kneaddataDB/human_genome_bowtie2 \
--metaphlan_db /home/mypool/projects/nasa_pipeline/metaphlanDB \
--humann_nucleotide_db /home/mypool/projects/nasa_pipeline/humannDB/chocophlan \
--humann_protein_db /home/mypool/projects/nasa_pipeline/humannDB/uniref \
--humann_pathway_db metacyc \
--goterm_db /home/mypool/projects/nasa_pipeline/gotermsDB \
--metaphlan_extra_analysis true \
--dev -resume
```

### Parameters

#### Global Options

| Parameter         | Default                | Description                                                                                  | Options                                  |
|-------------------|------------------------|----------------------------------------------------------------------------------------------|------------------------------------------|
| `-profile`        | `local`                | Execution profile and container engine                                                       | `local`, `slurm`, `conda`, `docker`, `singularity`, `apptainer` |
| `--samplesheet`   | NULL                   | CSV file with sample names and paired-end FASTQ file paths (R1, R2)                          | -                                        |
| `--output`        | ``${baseDir}/results`` | Output directory for results                                                                 | -                                        |
| `--dev`           | -                      | Run R code directly from the local GitHub branch instead of within the container (for development or debugging; recommended only if you are modifying or testing pipeline R scripts) | -                                        |

> [!NOTE]
> By default, when using container options (`docker`, `singularity`, or `apptainer`), all R code is executed within the respective container environment. If you wish to run the R code directly from the local GitHub branch (for development or debugging), add the `--dev` flag to your pipeline command. This is useful because some parts of the pipeline are still in development.

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
| `--goterm_time`             | -            | Process time (cluster only)                                      | -                              |
| `--cpa_cpu`                 | 2            | CPUs for consensus pathway analysis                              | -                              |
| `--cpa_mem`                 | 16 GB        | Memory for consensus pathway analysis                            | -                              |
| `--cpa_time`                | -            | Process time (cluster only)                                      | -                              |

#### Descriptive Profiling Analysis Parameters

| Parameter                   | Default      | Description                                                      | Options                        |
|-----------------------------|--------------|------------------------------------------------------------------|--------------------------------|
| `--humann_renorm_units`     | cpm          | Pathway abundance normalization units                            | `cpm`, `relab`                 |
| `--norm_path_cpu`           | 2            | CPUs for descriptive profiling                                   | -                              |
| `--norm_path_mem`           | 16 GB        | Memory for descriptive profiling                                 | -                              |
| `--norm_path_time`          | -            | Process time (cluster only)                                      | -                              |
| `--join_path_cpu`           | 1            | CPUs for join pathway                                            | -                              |
| `--join_path_mem`           | 8 GB         | Memory for join pathway                                          | -                              |
| `--join_path_time`          | -            | Process time (cluster only)                                      | -                              |
| `--desc_cpu`                | 1            | CPUs for descriptive analysis                                    | -                              |
| `--desc_mem`                | 8 GB         | Memory for descriptive analysis                                  | -                              |
| `--desc_time`               | -            | Process time (cluster only)                                      | -                              |

### Configuation Profiles

The pipeline currently supports two execution profiles set by `-profile` option by `,`. The options are:

- **local**
- **slurm**

And, four execution engines:

- **conda**
- **docker**
- **singularity**
- **apptainer**

> [!TIP]
> To run the pipeline in the `cluster` profile, please configure the `conf/cluster.config` file with your cluster settings.
> Accordingly, set the memory and time limits for each process in the `nextflow.config` file.


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

Please refer to the [wiki document](docs/wiki.md) for **Step 3**, **Step 4**, and **Step 5**. These steps guide you through downloading or providing the required **KneadData Database**, **MetaPhlAn Database**, and **HUMAnN Nucleotide and Protein Databases**.

## Output

The pipeline generates the following output files:

> [!NOTE]
> The output files are stored in nested directories within the user-defined output directory.

1. **kneaddata_out**
    -  Processed (decontaminated, trimmed) FASTQ files in `.fastq.gz` format. 
    -  `_kneaddata.log` file containing the log of the KNEADING_DATA process for each samples.
    -  `/fastqc` directory containing FastQC reports for the processed FASTQ files in `.html` and the corresponding `.zip` files.

> [!NOTE]
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
