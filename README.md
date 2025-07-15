# EPSCoR NASA Special Project: Pathway-level, consensus analysis of microbiome profiling data

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)


## Introduction

[Pipeline-name], is a scalable and reproducible pipeline designed for comprehensive taxonomic, functional and consensus pathway analysis of metagenomes. Implemented using Nextflow, the pipeline supports scalable execution across local, HPC, and cloud environments. It integrates established tools for quality control, taxonomic, and functional profiling, and introduces a novel consensus-based pathway enrichment strategy that combines multiple statistical methods to improve reliability across experimental conditions. When group comparisons are not defined, the pipeline offers a descriptive profiling mode that links microbial taxa with associated functions. We applied [Pipeline-name] to three distinct metagenomic datasets – including a space habitat study investigating samples from disrupted bioreactor samples, soil rhizosphere samples, and [description of Dr. Frese’s dataset] – demonstrating its versatility and ability to uncover biologically relevant functional patterns.

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

## Usage

> [!NOTE]
> If you are new to Nextflow, please refer to the [this page](https://nextflow.io/docs/latest/install.html) on how to set-up Nextflow.
> We have tested the pipeline on the above mentioned datasets. Please, download the datasets or choose your own dataset.

Create a samplesheet in CSV format. Below are two examples of how the samplesheet should look like:
  
  Example 1: For dataset with factor (condition) information.

  ```csv
  sample,factor,fastq_1,fastq_2
  SAMPLE1-ID,NN,sample1_R1.fastq.gz,sample1_R2.fastq.gz
  SAMPLE2-ID,NN,sample2_R1.fastq.gz,sample2_R2.fastq.gz
  SAMPLE3-ID,NN,sample3_R1.fastq.gz,sample3_R2.fastq.gz
  ```

  Example 2: For dataset without factor (condition) information.

  ```csv
  sample,fastq_1,fastq_2
  SAMPLE1-ID,sample1_R1.fastq.gz,sample1_R2.fastq.gz
  SAMPLE2-ID,sample2_R1.fastq.gz,sample2_R2.fastq.gz
  SAMPLE3-ID,sample3_R1.fastq.gz,sample3_R2.fastq.gz
  ```

1. Run the command below to execute the pipeline.

  ```bash
nextflow main.nf \\
            -profile <PROFILE> \\
            --samplesheet path/to/INPUT.csv \\
            --output /path/to/OUTPUT_DIR \\
            --kneaddata_db_path /path/to/folder/kneaddata/database/download/ \\
            --kneaddata_db REFERENCE_DB:BUILD \\
            --metaphlan_db_path /path/to/folder/metaphlan/database/download \\
            --humann_nucleotide_db NUCLEOTIDE_DB:BUILD \\
            --humann_nuc_db_path /path/to/folder/humann/nucleotide/databases/download \\
            --humann_protein_db <PROTEIN_DB:BUILD> \\
            --humann_prot_db_path </path/to/folder/humann/protein/databases/download> \\
            --humann_pathway_db <PATHWAY_DB>
```

### Parameter options

* `-profile` : Specify the profile to use for running the pipeline in local or on the HPC.
    This can be set to the following:
    1. _local_ ( To run the pipeline on local machine, uses separate conda environment for each process )
    2. _cluster_ ( To run the pipeline on HPC, uses separate conda container for each process )

* `--samplesheet` : Path to a CSV file where each row specifies the sample name and the file paths to paired-end FASTQ files (R1 and R2), separated by commas.

* `--output` : Path to the output directory where results will be saved.

* `--kneaddata_db_path` : Path to saved or where to save kneaddata databases (default: ./kneaddata_db/)

* `--kneaddata_db` : Comma seperated list with no spaces of databases for kneadata to use in database:build format  (default: human_transcriptome:bowtie2,ribosomal_RNA:bowtie2)

> [!Possible options for `--kneaddata_db`]
>
> - human_transcriptome:bowtie2
> - ribosomal_RNA:bowtie2
> - mouse_C57BL:bowtie2
> - dog_genome:bowtie2
> - cat_genome:bowtie2
> - human_genome:bmtagger


> [!IMPORTANT]
> There currently is a bug in the current *kneaddata v0.12.2* release of kneaddata for human_genome:bowtie2. This will be fixed in _0.12.3_. You can manually create the *"human_genome_bowtie2"* directory and manually download the correct file here `https://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_hg39_T2T_Bowtie2_v0.1.tar.gz` and extract it into the *kneaddata_path/human_genome_bowtie2* directory.
> 
> Then, within that directory run the following command: `touch .done`

* `--metaphlan_db_path` : Path to the directory where metaphlan databases are saved or will be downloaded (default: ./metaphlan_db/)

* `--humann_nucleotide_db` : Comma separated list with no spaces of nucleotide databases for humann to use in database:build format (default: chocophlan:full)
    Possible options:
        - chocophlan:full

* `--humann_nuc_db_path` : Path to the directory where humann nucleotide database is saved or will be downloaded (default: ./humann_nucdb/)

* `--humann_protein_db` : Comma separated list with no spaces of protein databases for humann to use in database:build format (default: uniref:uniref50_diamond)

> [!Possible options for `--humann_protein_db`
> 
> - uniref:uniref50_diamond
> - uniref:uniref90_diamond
> - uniref:uniref50_ec_filtered_diamond
> - uniref:uniref90_ec_filtered_diamond

* `--humann_prot_db_path` : Path to the directory where humann protein database is saved or will be downloaded (default: ./humann_protdb/)

* `--humann_pathway_db` : Specify the database to use for pathway {metacyc, unipathways} computations (default: metacyc)

## Reporting

## Pipeline Output

## Credits

## Contribution and Support

## Citations
