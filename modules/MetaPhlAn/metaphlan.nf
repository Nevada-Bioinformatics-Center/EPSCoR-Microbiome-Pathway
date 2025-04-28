/*
 * Module: MetaPhlAn
 * Description: This module performs MetaPhlAn analysis on metagenomic data.
 * Version: v4.1.1
 * Author: Kanishka Manna
 * Date: 2025-04-11
 */


process TAXONOMIC_PROFILING {
    tag "MetaPhlAn"
    publishDir "${params.output}/metaphlan_out", mode: 'copy'

    input:
        tuple val(sample_id), path(reads)

    output:
        path( "${sample_id}_profile.tsv" ),  emit: profile_taxa
        path( "${sample_id}_metaphlan.log"), emit: metaphlan_log

    script:
    """
    # Download MetaPhlAn database if not already present
    # metaphlan --install --bowtie2db ${params.metaphlan_db}

    # Run MetaPhlAn on the input reads
    metaphlan \
    ${reads} \
    --input_type fastq \
    --output_file ${sample_id}_profile.tsv \
    --bowtie2out ${sample_id}_bowtie2.bz2 \
    --nproc 4 \
    &> ${sample_id}_metaphlan.log
    """
}