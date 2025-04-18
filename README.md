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
            --samplesheet INPUT.csv \\
            --output OUTPUT_DIR \\
            --library PE \\
            --kneaddata_db REFERENCE_DB \\
            --decontaminate_pairs {strict, lenient, unpaired}
```

## Parameter options

`--samplesheet` : path to a CSV file where each row specifies the sample name and the file paths to paired-end FASTQ files (R1 and R2), separated by commas.

`--output` : path to the output directory where results will be saved.

`--library` : sequencing library type, e.g., PE for paired-end or SE for single-end. (default: PE)

`--kneaddata_db` : path to the KneadData reference database required for decontamination.

`--decontaminate_pairs` : specify the decontamination mode, either 'strict', 'lenient' or 'unpaired'. (default: lenient)

## Note

The pipeline is currently running till, QC step.

## Things to do

1. Fix KneadData module - currently, KneadData is not able to call Bowtie2
2. Add a database/reference module for KneadData, MetaPhlAn and HUMAnN
3. Edit MetaPhlAn
4. Edit HUMAnN
5. Add Pathway analysis
6. Add Visualization
