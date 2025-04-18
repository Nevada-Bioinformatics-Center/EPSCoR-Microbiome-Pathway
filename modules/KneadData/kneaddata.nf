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
        tuple val(sample_id), path(read)
        path(params.kneaddata_db)

    output:
        path( "*.fastq.gz" ),   emit: kneaddata_fastqc
        path( "*.log" ),        emit: kneaddata_log

    script:
        """
        kneaddata \
        --input1 ${read[0]} \
        --input2 ${read[1]} \
        --reference-db ${params.kneaddata_db} \
        --output  "kneaddata_out" \
        --decontaminate-pairs ${params.decontaminate_pairs} \
        --log ${sample_id}_kneaddata.log
        """
}