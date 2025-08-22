/*
 * Module: MetaPhlAn
 * Description: This module performs MetaPhlAn analysis on metagenomic data.
 * Version: v4.1.1
 * Author: Kanishka Manna and Hans Vasquez-Gross
 */

process TAXONOMIC_PROFILING {

    tag "MetaPhlAn: ${sample_id}"

    label 'metaphlan'

    publishDir "${params.output}/metaphlan_out", mode: 'copy'

    input:
        tuple val(sample_id), path(merged_fastq), path(metaphlan_db_dir)

    output:
        path( "${sample_id}_profile.tsv" ), emit: profiled_taxa
        path( "${sample_id}_profile.txt" ), emit: profiled_taxa_txt
        path( "${sample_id}_metaphlan.log" ), emit: metaphlan_log

    script:
    """
    metaphlan \\
        ${merged_fastq} \\
        --input_type fastq \\
        --bowtie2db ${metaphlan_db_dir} \\
        --output_file ${sample_id}_profile.tsv \\
        --bowtie2out ${sample_id}_bowtie2.bz2 \\
        --unclassified_estimation \\
        -t rel_ab_w_read_stats \\
        --nreads 100 \\
        --nproc ${params.metaphlan_nproc} \\
        --offline \\
        --index ${params.metaphlan_index} \\
        ${params.metaphlan_extra} \\
        &> ${sample_id}_metaphlan.log
    
    cp ${sample_id}_profile.tsv ${sample_id}_profile.txt
    """
}
