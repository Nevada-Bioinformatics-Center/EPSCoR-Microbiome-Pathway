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
        path(params.metaphlan_db)

    output:
        path( "${sample_id}_taxa-profile.tsv" ),  emit: metaphlan_profile

    script:
    """
    metaphlan $reads \\
        --input_type fastq \\
        --bowtie2db ${params.metaphlan_db} \\
        --output_file ${sample_id}_profile.tsv
    """
}