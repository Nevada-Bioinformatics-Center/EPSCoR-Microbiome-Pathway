/*
 * Module: HUMAnN
 * Description: This module downloads the nucleotide & protein databases for HUMAnN analysis.
 * Version: v3.9
 * Author: Kanishka Manna
 * Date: 2025-05-08
 */

process DOWNLOAD_HUMANN_NUCLEOTIDE_DB {

    tag "humann_database_download: ${nuc_database};${nuc_build}"

    label 'humann_conda'
    label 'medium'

    input:
    tuple val(nuc_database), val(nuc_build), val(install_dir)

    output:
    val("${install_dir}/${nuc_database}"), emit: downloaded_nucleotide_db  // Output just the path string

    script:
    """
    if [ ! -f "${install_dir}/.done" ]; then
        echo "Downloading ${nuc_database} ${nuc_build} to ${install_dir}"
        mkdir -p ${install_dir}
        humann_databases --download ${nuc_database} ${nuc_build} ${install_dir}
        touch ${install_dir}/.done
    else
        echo "Database already exists at ${install_dir}, skipping download."
    fi
    """
}


process DOWNLOAD_HUMANN_PROTEIN_DB {

    tag "humann_database_download: ${prot_database};${prot_build}"

    label 'humann_conda'
    label 'medium'

    input:
    tuple val(prot_database), val(prot_build), val(install_dir)

    output:
    val("${install_dir}/${prot_database}"), emit: downloaded_nucleotide_db  // Output just the path string

    script:
    """
    if [ ! -f "${install_dir}/.done" ]; then
        echo "Downloading ${prot_database} ${prot_build} to ${install_dir}"
        mkdir -p ${install_dir}
        humann_databases --download ${prot_database} ${prot_build} ${install_dir}
        touch ${install_dir}/.done
    else
        echo "Database already exists at ${install_dir}, skipping download."
    fi
    """
}