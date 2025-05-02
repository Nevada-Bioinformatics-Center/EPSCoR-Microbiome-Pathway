/*
 * Module: MetaPhlAn
 * Description: This module performs MetaPhlAn analysis on metagenomic data.
 * Version: v4.1.1
 * Author: Kanishka Manna and Hans Vasquez-Gross
 * Date: 2025-04-11
 */

process DOWNLOAD_METPHLAN_DB {
    tag "metaphlan_db_download"

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

process TAXONOMIC_PROFILING {
    tag { "MetaPhlAn: ${sample_id}" }  // ✅ correct way to defer interpolation
    publishDir "${params.output}/metaphlan_out", mode: 'copy'


    input:
        tuple val(sample_id), path(reads), path(metaphlan_db_dir)

    output:
        path("${sample_id}_profile.tsv"), emit: profiled_taxa
        path("${sample_id}_metaphlan.log"), emit: metaphlan_log

    script:
    """
    export DEFAULT_DB_FOLDER=${metaphlan_db_dir}

    metaphlan \\
        ${reads[0]} ${reads[1]} \\
        --input_type fastq \\
        --bowtie2db \$DEFAULT_DB_FOLDER \\
        --output_file ${sample_id}_profile.tsv \\
        --bowtie2out ${sample_id}_bowtie2.bz2 \\
        --nproc 4 \\
        &> ${sample_id}_metaphlan.log
    """
}
