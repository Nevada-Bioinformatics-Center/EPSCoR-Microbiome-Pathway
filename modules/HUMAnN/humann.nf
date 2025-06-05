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
    label 'humann_docker'
    label 'high'
    
    publishDir "${params.output}/humann_out", mode: 'copy'

    input:
        tuple val(sample_id), path(reads), val(nucleotide_db), val(protein_db), path(taxonomic_profile)

    output:
        path( "${sample_id}_genefamilies.tsv" ),  emit: gene_fam
        path( "${sample_id}_pathabundance.tsv" ), emit: path_abundance
        path( "${sample_id}_pathcoverage.tsv" ),  emit: path_coverage
        path( "${sample_id}_humann.log" ),    emit: humann_log

    script:
        def nucleotide_db_flags = nucleotide_db.collect { "--nucleotide-database ${it}" }.join(' ')
        
        def protein_db_flags = protein_db.collect { "--protein-database ${it}" }.join(' ')
    """
    # Run HUMAnN on the input reads
    humann \\
        --input ${reads} \\
        --input-format fastq \\
        ${nucleotide_db_flags} \\
        ${protein_db_flags} \\
        --taxonomic-profile ${taxonomic_profile} \\
        --pathways ${params.humann_pathway_db} \\
        --output "./" \\
        --output-format tsv \\
        --remove-temp-output \\
        --threads 5 \\
        &> ${sample_id}_humann.log
    """
}