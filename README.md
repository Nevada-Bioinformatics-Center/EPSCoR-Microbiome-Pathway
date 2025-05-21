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
            --humann_db_path /path/to/folder/humann/databases/download \\
            --humann_nucleotide_db NUCLEOTIDE_DB:BUILD \\
            --humann_protein_db PROTEIN_DB:BUILD \\
            --humann_pathway_db PATHWAY_DB
```

## Parameter options
