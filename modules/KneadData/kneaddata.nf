/*
 * Module: KneadData
 * Description: This module performs kneading of metagenomic data.
 * Version: v0.12.2
 * Author: Kanishka Manna
 * Date: 2025-04-09
 */

process KNEADING_DATA {
    tag "KneadData"
    publishDir "${params.output}/kneaddata_out", mode: 'copy'

    input:
        tuple val(sample_id), path(read), val(kneaddata_db_paths) 

    output:
        path( "*.fastq" ),   emit: kneaddata_fastq
        path( "*.log" ),        emit: kneaddata_log
        path( "fastqc/${sample_id}*.html" ),        emit: kneaddata_fastqc_html
        path( "fastqc/${sample_id}*.zip" ),        emit: kneaddata_fastqc_zip



    script:
        def db_flags = kneaddata_db_paths.collect { "--reference-db ${it}" }.join(' ')
        """
        kneaddata \
        --input1 ${read[0]} \
        --input2 ${read[1]} \
        --threads 10 \
        --processes 2 \
        ${db_flags} \
        --output  "./" \
        --output-prefix ${sample_id} \
        --decontaminate-pairs ${params.decontaminate_pairs} \
        --run-fastqc-start \
        --run-fastqc-end \
        --log ${sample_id}_kneaddata.log
        """
}