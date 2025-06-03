# EPSCoR NASA Special Project: Pathway-level, consensus analysis of microbiome profiling data

## Download the data

Download the data from [**Open Science Data Repository (OSDR)**](https://www.nasa.gov/osdr/). We are currently using the [_OSD-286_](https://osdr.nasa.gov/bio/repo/data/studies/OSD-286) dataset.
Description: Metagenome analysis of ISS cargo resupply vehicles (CRV).
Factors: 'Spaceflight', 'Treatment', 'Sample Location'
Assay: Metagenomic sequencing - Whole Genome Shotgun Sequencing - Illumina Nextera Kit
Device Platform: Illumina
Samples: 8
Sample Type: Paired-end

## To run the pipeline

- Install Nextflow on your system.
- Create a conda environment and install the **"Biobakery"** packages - **"KneadData"**, **"MetaPhlAn"**, **"HUMAnN"** and **"FastQC"**.
- Use VSCode for better nextflow editing experience.

```
nextflow main.nf \\
            -profile local \\
            --samplesheet path/to/INPUT.csv \\
            --output /path/to/OUTPUT_DIR \\
            --kneaddata_db_path /path/to/folder/kneaddata/database/download/ \\
            --kneaddata_db REFERENCE_DB:BUILD \\
            --metaphlan_db_path /path/to/folder/metaphlan/database/download \\
            --humann_nucleotide_db NUCLEOTIDE_DB:BUILD \\
            --humann_nuc_db_path /path/to/folder/humann/nucleotide/databases/download \\
            --humann_protein_db PROTEIN_DB:BUILD \\
            --humann_prot_db_path /path/to/folder/humann/protein/databases/download \\
            --humann_pathway_db PATHWAY_DB
```

## Parameter options

* `--samplesheet` : Path to a CSV file where each row specifies the sample name and the file paths to paired-end FASTQ files (R1 and R2), separated by commas.

* `--output` : Path to the output directory where results will be saved.

* `--kneaddata_db_path` : Path to saved or where to save kneaddata databases (default: ./kneaddata_db/)

* `--kneaddata_db` : Comma seperated list with no spaces of databases for kneadata to use in database:build format  (default: human_transcriptome:bowtie2,ribosomal_RNA:bowtie2)

    Possible options:
        - human_transcriptome:bowtie2
        - ribosomal_RNA:bowtie2
        - mouse_C57BL:bowtie2
        - dog_genome:bowtie2
        - cat_genome:bowtie2
        - human_genome:bmtagger

    There currently is a bug in the current *kneaddata v0.12.2* release of kneaddata for human_genome:bowtie2. This will be fixed in _0.12.3_. You can manually create the *"human_genome_bowtie2"* directory and manually download the correct file here 
    `https://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_hg39_T2T_Bowtie2_v0.1.tar.gz` and extract it into the *kneaddata_path/human_genome_bowtie2* directory. The within that directory run the following command: touch `.done`

* `--metaphlan_db_path` : Path to the directory where metaphlan databases are saved or will be downloaded (default: ./metaphlan_db/)

* `--humann_nucleotide_db` : Comma separated list with no spaces of nucleotide databases for humann to use in database:build format (default: chocophlan:DEMO)
    Possible options:
        - chocophlan:full
        - chocophlan:DEMO

* `--humann_nuc_db_path` : Path to the directory where humann nucleotide database is saved or will be downloaded (default: ./humann_nucdb/)

* `--humann_protein_db` : Comma separated list with no spaces of protein databases for humann to use in database:build format (default: uniref:DEMO_diamond)

    Possible options:
        - uniref:uniref50_diamond
        - uniref:uniref90_diamond
        - uniref:uniref50_ec_filtered_diamond
        - uniref:uniref90_ec_filtered_diamond
        - uniref:DEMO_diamond

* `--humann_prot_db_path` : Path to the directory where humann protein database is saved or will be downloaded (default: ./humann_protdb/)

* `--humann_pathway_db` : Specify the database to use for pathway {metacyc, unipathways} computations (default: metacyc)
