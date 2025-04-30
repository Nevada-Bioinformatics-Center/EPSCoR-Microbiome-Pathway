/*
 * Module: KneadData
 * Description: This module performs kneading of metagenomic data.
 * Version: v0.12.2
 * Author: Hans Vasquez-Gross
 * Date: 2025-04-18
 */


process download_kneaddata_db {

    tag "kneaddata_download: ${database};${build}"

    input:
    tuple val(database), val(build), val(install_dir)

    output:
    val(install_dir), emit: downloaded_db  // Output just the path string

    script:
    """
    if [ ! -f "${install_dir}/.done" ]; then
        echo "Downloading ${database} ${build} to ${install_dir}"
        mkdir -p ${install_dir}
        kneaddata_database --download ${database} ${build} ${install_dir}
        touch ${install_dir}/.done
    else
        echo "Database already exists at ${install_dir}, skipping download."
    fi
    """
}