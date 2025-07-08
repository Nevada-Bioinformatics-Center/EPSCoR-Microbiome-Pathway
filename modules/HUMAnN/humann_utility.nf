/*
 * Module: HUMAnN Utility tools
 * Description: This module performs HUMAnN analysis on metagenomic data.
 * Version: v3.9
 * Author: Kanishka Manna
 * Date: 2025-04-28
 */


//  This process will join (merge) multiple single-sample output files 
//  into a single table with multi-samples
<<<<<<< Updated upstream
<<<<<<< Updated upstream
//process join_tables {
//    tag "humann_join_tables"
//    publishDir "${params.output}/humann_out", mode: 'copy'

//    input:
//        path(genefamilies)
//        path(pathabundance)
//        path(pathcoverage)


//    output:
//        path( "joined_genefamilies.tsv" ),  emit: merged_genefamilies
//        path( "joined_pathabundance.tsv" ), emit: merged_pathabundance
//        path( "joined_pathcoverage.tsv" ), emit: merged_pathcoverage

//    script:
//    """
//    humann_join_tables --input ${genefamilies} --output ${merged_genefamilies}
//    humann_join_tables --input ${pathabundance} --output ${merged_pathabundance}
//    humann_join_tables --input ${pathcoverage} --output ${merged_pathcoverage}
//    """
//}
=======
=======
>>>>>>> Stashed changes
process NORMALIZE_PATHWAY_ABUNDANCE {
    tag "humann_renorm_table"

    label 'humann_conda'

    publishDir "${params.output}/humann_out/pathabundance", mode: 'copy'

    input:
        tuple val (sample_id), path(pathabundance)

    output:
        path( "${sample_id}_renorm_pathabundance.tsv" ), emit: renorm_pathabundance
    script:
    """
    humann_renorm_table --input ${pathabundance} \
    --output "${sample_id}_renorm_pathabundance.tsv" \
    --units relab
    """
}
>>>>>>> Stashed changes


// This will normalize the merged pathway abundance to relative abundances
process JOIN_PATHWAY_ABUNDANCE {
    tag "humann_join_tables"
<<<<<<< Updated upstream

    label 'humann_conda'

    publishDir "${params.output}/humann_out/merge", mode: 'copy'

    input:
        path(renorm_pathabundance_dir)

    output:
        path( "joined_norm_pathabundance.tsv" ), emit: join_table

=======

    label 'humann_conda'

    publishDir "${params.output}/humann_out/merge", mode: 'copy'

    input:
        path(renorm_pathabundance_dir)

    output:
        path( "joined_norm_pathabundance.tsv" ), emit: join_table

>>>>>>> Stashed changes
    script:
    """
    humann_join_tables --input ${renorm_pathabundance_dir} --output joined_norm_pathabundance.tsv
    """
}



// This process will generate descriptive statistics for the pathway abundance
process DESC_PROFILING {
    tag "Descriptive Profiling"

    label 'desc_conda'

    publishDir "${params.output}/humann_out/desc", mode: 'copy'

    input:
        path(joined_pathabundance)
    
    output:
        path("top_genus_out.tsv"), emit: top_genus
        path("top_species_out.tsv"), emit: top_species
    
    script:
    """
    desc_prof.R ${joined_pathabundance} top_genus_out.tsv top_species_out.tsv
    """
}