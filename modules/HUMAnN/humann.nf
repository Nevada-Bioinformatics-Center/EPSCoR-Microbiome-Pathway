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
        path(params.humann_nucleotide_db)
        path(params.humann_protein_db)
        path(profiled_taxa)


    output:
        path( "${sample_id}/*_genefamilies.tsv" ),  emit: humann_genefamilies
        path( "${sample_id}/*_pathabundance.tsv" ), emit: humann_pathabundance
        path( "${sample_id}/*_pathcoverage.tsv" ),  emit: humann_pathcoverage
        path( "${sample_id}/*_humann.log" ),    emit: humann_log

    script:
    """
    humann \\
        --input ${reads} \\
        --input-format fastq \\
        --nucleotide-database ${params.humann_nucleotide_db} \\
        --protein-database ${params.humann_protein_db} \\
        --taxonomic-profile ${profiled_taxa} \\
        --output ${sample_id}_humann_out \\
        --output-basename ${sample_id} \\
        --output-format tsv \\
        --remove-temp-output
    """
}