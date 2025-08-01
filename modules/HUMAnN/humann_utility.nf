/*
 * Module: HUMAnN Utility tools
 * Description: This module performs HUMAnN analysis on metagenomic data.
 * Version: v3.9
 * Author: Kanishka Manna
 */


// This will normalize the pathway abundance to relative abundances
process NORMALIZE_PATHWAY_ABUNDANCE {
    tag "humann_renorm_table"

    label 'humann_conda'

    publishDir "${params.output}/humann_out/pathabundance", mode: 'symlink'

    input:
        tuple val (sample_id), path(pathabundance)

    output:
        path( "${sample_id}_renorm_pathabundance.tsv" ), emit: renorm_pathabundance
    script:
    """
    humann_renorm_table --input ${pathabundance} \
    --output "${sample_id}_renorm_pathabundance.tsv" \
    --units ${params.humann_renorm_units}
    """
}



// This will normalize the merged pathway abundance to relative abundances
process JOIN_PATHWAY_ABUNDANCE {
    tag "humann_join_tables"

    label 'humann_conda'

    publishDir "${params.output}/humann_out/merge", mode: 'copy'

    input:
        path(renorm_pathabundance_dir)

    output:
        path( "joined_norm_pathabundance.tsv" ), emit: join_table

    script:
    """
    humann_join_tables --input ${renorm_pathabundance_dir} --output joined_norm_pathabundance.tsv
    """
}



// This process will generate descriptive statistics for the pathway abundance
process DESCRIPTIVE_PROFILING {
    tag "Descriptive Profiling"

    label 'desc_conda'

    publishDir "${params.output}/desc_out", mode: 'copy'

    input:
        path(joined_pathabundance)
    
    output:
        path("top_path_taxa_results.tsv"), emit: top_taxa
    
    script:
    """
    desc_prof.R ${joined_pathabundance} top_path_taxa_results.tsv
    """
}

