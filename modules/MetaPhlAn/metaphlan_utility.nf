/*
 * Module: MetaPhlAn
 * Description: This module performs MetaPhlAn analysis on metagenomic data.
 * Version: v4.1.1
 * Author: Kanishka Manna and Hans Vasquez-Gross
 */

// Process to run MetaPhlAn merging of taxonomic profiling tables
process MERGE_TAXONOMIC_TABLES {
    tag "Merge MetaPhlAn tables"
    
    label 'metaphlan_conda'

    publishDir "${params.output}/metaphlan_out/merge", mode: 'symlink'

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

    publishDir "${params.output}/metaphlan_out/plots", mode: 'move'

    input:
        path(merged_taxa_profile) 
        path(samplesheet_file)

    output:
        path("*.png"), emit: taxa_viz

    script:
    """
    taxa_viz.r ${merged_taxa_profile} ${samplesheet_file} ${params.output}/metaphlan_out/plots
    """
}