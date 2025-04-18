/*
 * Module: HUMAnN
 * Description: This module performs HUMAnN analysis on metagenomic data.
 * Version: v3.9
 * Author: Kanishka Manna
 * Date: 2025-04-09
 */

process FUNCTIONAL_PROFILING {
    tag "HUMAnN"
    publishDir "${params.output}/humann_out", mode: 'copy'

    input:
        tuple val(sample_id), path(reads)
        path(params.humann_db)

    output:
        path( "${sample_id}/*_genefamilies.tsv" ),  emit: humann_genefamilies
        path( "${sample_id}/*_pathabundance.tsv" ), emit: humann_pathabundance
        path( "${sample_id}/*_pathcoverage.tsv" ),  emit: humann_pathcoverage
        path( "${sample_id}/*_HUMANnN.log" ),       emit: humann_log

    script:
    """
    humann \\
        --input ${reads} \\
        --output ${sample_id}_humann_out \\
        --output-prefix ${sample_id} \\
        --database ${params.humann_db} \\
        --sequencer-source ${params.sequencer_source} \\
        --decontaminate-pairs ${params.decontaminate_pairs} \\
        --remove-intermediate-output \\
        --reorder
    """
}