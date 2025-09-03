/*
 * Module: MetaPhlAn Utility tool and Visualization
 * Description: This module performs downstream analysis on merged taxonomic file from MetaPhlAn.
 * Version: v4.1.1
 * Author: Kanishka Manna
 */

// Process to run MetaPhlAn merging of taxonomic profiling tables
process MERGE_TAXONOMIC_TABLES {
    tag "Merge MetaPhlAn tables"
    
    label 'metaphlan_conda'

    publishDir "${params.output}/metaphlan_out/merge", mode: 'copy'

    input:
        path(metaphlan_tables)
    
    output:
        path("merged_taxa_profile.tsv"), emit: merged_taxa_profile

    script:
    """
    merge_metaphlan_tables.py ${metaphlan_tables} > merged_taxa_profile.tsv
    """
}


// Process to run Taxonomic visualization analysis
process TAXONOMIC_VISUALIZATION {
    tag "MetaPhlAn Visualization"

    label 'taxviz_conda'

    publishDir "${params.output}/metaphlan_out", mode: 'move'

    input:
        path(merged_taxa_profile) 
        path(samplesheet_file)

    output:
        path("*"), emit: taxa_viz

    script:
    """
    taxa_viz.r ${merged_taxa_profile} ${samplesheet_file}
    """
}