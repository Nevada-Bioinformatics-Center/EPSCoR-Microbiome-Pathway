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
            --samplesheet path/to/INPUT.csv \\
            --output /path/to/OUTPUT_DIR \\
            --library PE \\
            --kneaddata_path /path/to/folder/kneaddata/database/download/ \\
            --kneaddata_db REFERENCE_DB \\
            --nucleotide_db /path/to/downloaded/nucleotide/database \\
            --protein_db /path/to/downloaded/protein/database \\
            --pathway_db metacyc
```

## Parameter options

* `--samplesheet` : Absolute path to a CSV file where each row specifies the sample name and the file paths to paired-end FASTQ files (R1 and R2), separated by commas.

* `--output` : Absoulute path to the output directory where results will be saved.

* `--library` : sequencing library type, e.g., PE for paired-end or SE for single-end. (default: PE)

* `--kneaddata_path` : Absolute path to saved or where to save knead_data databases (default: ./kneaddata_db/)

* `--kneaddata_db` : comma seperated list with no spaces of databases for kneadata to use in database:build format  (default: human_transcriptome:bowtie2,ribosomal_RNA:bowtie2) 
    Possible options: 
        `human_transcriptome:bowtie2`, `ribosomal_RNA:bowtie2`, `mouse_C57BL:bowtie2`, `dog_genome:bowtie2`, `cat_genome:bowtie2`, `human_genome:bmtagger`

        There currently is a bug in the current kneaddata v0.12.2 release of kneaddata for human_genome:bowtie2. This will be fixed in 0.12.3.
        You can manually create the "human_genome_bowtie2" directory and manually download the correct file here 
        https://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_hg39_T2T_Bowtie2_v0.1.tar.gz

        and extract it into the kneaddata_path/human_genome_bowtie2 directory. The within that directory run the following command: touch .done
            
* `--nucleotide_db` : Specify the directory containing the nucleotide database

* `--protein_db` : Specify the directory containing the protein database
            
* `--pathway_db` : Specify the database to use for pathway {metacyc, unipathways} computations (default: metacyc)