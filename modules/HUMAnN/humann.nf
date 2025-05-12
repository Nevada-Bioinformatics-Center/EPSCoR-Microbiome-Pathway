/*
 * Module: HUMAnN
 * Description: This module performs HUMAnN analysis on metagenomic data.
 * Version: v3.9
 * Author: Kanishka Manna
 * Date: 2025-04-09
 */

process FUNCTIONAL_PROFILING {
    tag "HUMAnN: ${sample_id}"

    label 'humann_conda'
    label 'high'
    
    publishDir "${params.output}/humann_out", mode: 'copy'

    input:
        tuple val(sample_id), path(reads), path(taxonomic_profile)
        path(nucleotide_db)
        path(protein_db)

    output:
        path( "${sample_id}_genefamilies.tsv" ),  emit: gene_fam
        path( "${sample_id}_pathabundance.tsv" ), emit: path_abundance
        path( "${sample_id}_pathcoverage.tsv" ),  emit: path_coverage
        path( "${sample_id}_humann.log" ),    emit: humann_log

    script:
    """
    # Run HUMAnN on the input reads
    humann \\
        --input ${reads} \\
        --input-format fastq \\
        --nucleotide-database ${nucleotide_db} \\
        --protein-database ${protein_db} \\
        --taxonomic-profile ${taxonomic_profile} \\
        --pathways ${params.pathway_db} \\
        --output "./" \\
        --output-format tsv \\
        --remove-temp-output \\
        --threads 5 \\
        &> ${sample_id}_humann.log
    """
}