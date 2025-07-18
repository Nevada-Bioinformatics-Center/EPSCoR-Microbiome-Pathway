/*
 * Module: KneadData
 * Description: This module dowload database for kneading of metagenomic data.
 * Version: v0.12.2
 * Author: Hans Vasquez-Gross & Kanishka Manna
 */


process DOWNLOAD_KNEADDATA_DB {

    tag "kneaddata_download: ${database};${build}"

    label 'kneaddata_conda'
    label 'kneaddata_docker'
    label 'medium'

    publishDir "${params.database}/kneaddata_db", mode: 'symlink'

    input:
        tuple val(database), val(build)

    output:
        val("${params.database}/kneaddata_db/${database}_${build}"), emit: downloaded_db  // Output just the path string

    when:
        ! file("${params.database}/kneaddata_db/${database}_${build}").exists()

    script:
    """
    db_dir="${database}_${build}"

    mkdir -p "\$db_dir"
    
    kneaddata_database --download ${database} ${build} "\$db_dir"
    """
}