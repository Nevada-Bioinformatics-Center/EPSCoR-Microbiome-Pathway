/*
 * Module: HUMAnN
 * Description: This module downloads the nucleotide & protein databases for HUMAnN analysis.
 * Version: v3.9
 * Author: Kanishka Manna & Hans Vasquez-Gross
 * Date: 2025-05-08
 */

process DOWNLOAD_HUMANN_NUCLEOTIDE_DB {

    tag "humann_database_download: ${nuc_database};${nuc_build}"

    label 'humann'
    label 'medium'

    publishDir "${params.output}/humann_nucdb", mode: 'symlink'

    input:
        tuple val(nuc_database), val(nuc_build)

    output:
        path("${nuc_database}_${nuc_build}"), emit: downloaded_nucleotide_db  // Output just the path string

    when:
        ! file("${params.database}/humann_nucdb/${nuc_database}_${nuc_build}").exists()

    script:
    """
    db_dir="${nuc_database}_${nuc_build}"

    mkdir -p "\$db_dir"

    humann_databases --download ${nuc_database} ${nuc_build} "\$db_dir"
    """
}


process DOWNLOAD_HUMANN_PROTEIN_DB {

    tag "humann_database_download: ${prot_database};${prot_build}"

    label 'humann'
    label 'medium'

    publishDir "${params.output}/humann_protdb", mode: 'copy'

    input:
        tuple val(prot_database), val(prot_build)

    output:
        val("${prot_database}_${prot_build}"), emit: downloaded_nucleotide_db  // Output just the path string

    when:
        ! file("${params.database}/humann_protdb/${prot_database}_${prot_build}").exists()

    script:
    """
    db_dir="${prot_database}_${prot_build}"

    mkdir -p "\$db_dir"
        
    humann_databases --download ${prot_database} ${prot_build} "\$db_dir"
    """
}
