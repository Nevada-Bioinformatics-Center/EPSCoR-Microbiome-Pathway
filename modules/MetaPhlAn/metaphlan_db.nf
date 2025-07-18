/*
 * Module: MetaPhlAn
 * Description: This module performs MetaPhlAn analysis on metagenomic data.
 * Version: v4.1.1
 * Author: Hans Vasquez-Gross & Kanishka Manna
 */

process DOWNLOAD_METAPHLAN_DB {
    tag "metaphlan_db_download"
    
    label 'metaphlan_conda'
    label 'metaphlan_docker'
    label 'medium'

    publishDir "${params.database}", mode: 'symlink'

    output:
        val("${params.database}/metaphlan_db"), emit: metaphlan_db_dir

    when:
        ! file("${params.database}/metaphlan_db").exists()

    script:
    """
    db_dir="metaphlan_db"

    mkdir -p \$db_dir

    metaphlan --install --db_dir "\$db_dir" &> "\$db_dir/metaphlan_db.log"
    """
}