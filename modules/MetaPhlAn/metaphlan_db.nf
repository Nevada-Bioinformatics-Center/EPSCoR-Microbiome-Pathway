/*
 * Module: MetaPhlAn
 * Description: This module performs MetaPhlAn analysis on metagenomic data.
 * Version: v4.1.1
 * Author: Hans Vasquez-Gross
 * Date: 2025-04-11
 */

process DOWNLOAD_METPHLAN_DB {
    tag "metaphlan_db_download"
    label 'metaphlan_conda'

    input:
        val(install_dir)


    output:
        val(install_dir), emit: metaphlan_db_dir

    script:
    """
    if [ ! -f ${install_dir}/.done ]; then
        echo Downloading metaphlan database to ${install_dir}
        mkdir -p ${install_dir}
        export DEFAULT_DB_FOLDER=${install_dir}

        metaphlan --install --force_download --bowtie2db ${install_dir} &> ${install_dir}/metaphlan_db.log
        touch ${install_dir}/.done
    else
        echo "Database already exists at ${install_dir}, skipping download."
    fi
    
    """
}