/*
 * Module: MetaPhlAn
 * Description: This module performs MetaPhlAn analysis on metagenomic data.
 * Version: v4.1.1
 * Author: Kanishka Manna and Hans Vasquez-Gross
 * Date: 2025-04-11
 */

process TAXONOMIC_PROFILING {

    tag "MetaPhlAn: ${sample_id}"

    label 'metaphlan_conda'
    label 'metaphlan_docker'
    label 'high'

    publishDir "${params.output}/metaphlan_out", mode: 'copy'

    input:
        tuple val(sample_id), path(merged_fastq), path(metaphlan_db_dir)

    output:
        path("${sample_id}_profile.tsv"), emit: profiled_taxa
        path("${sample_id}_metaphlan.log"), emit: metaphlan_log

    script:
    """
    export DEFAULT_DB_FOLDER=${metaphlan_db_dir}

    metaphlan \\
        ${merged_fastq} \\
        --input_type fastq \\
        --bowtie2db \$DEFAULT_DB_FOLDER \\
        --output_file ${sample_id}_profile.tsv \\
        --bowtie2out ${sample_id}_bowtie2.bz2 \\
        --nproc 2 \\
        &> ${sample_id}_metaphlan.log
    """
}
