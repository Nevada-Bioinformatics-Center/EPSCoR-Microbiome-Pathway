/*
 * Module: HUMAnN Utility tools
 * Description: This module performs HUMAnN analysis on metagenomic data.
 * Version: v3.9
 * Author: Kanishka Manna
 * Date: 2025-04-28
 */


//  This process will join (merge) multiple single-sample output files 
//  into a single table with multi-samples
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


// This will normalize the merged pathway abundance to relative abundances
//process renorm_joined_table {
//    tag "humann_renorm_table"
//    publishDir "${params.output}/humann_out", mode: 'copy'

//    input:
//        path(merged_pathabundance)

//    output:
//        path( "joined_norm_pathabundance.tsv" ), emit: norm_table

//    script:
//    """
//    humann_renorm_table --input ${merged_pathabundance} --units relab --output ${norm_table}
//    """
//}