/*
 * Module: FastQC
 * Description: This module runs FastQC on paired-end FASTQ files.
 * Version: v0.12.1
 * Author: Kanishka Manna
 * Date: 2025-04-09
 */

process QUALITY_CONTROL {
    tag "QC"
    publishDir "${params.output}/fastqc_out/", mode: 'copy'

    input:
        tuple val(sample_id), path(reads)

    output:
        path( "${sample_id}*.html" ) , emit: fastqc_html
        path( "${sample_id}*_fastqc.zip" ) , emit: fastqc_zip


    script:
        """
        fastqc ${reads}
        """
}