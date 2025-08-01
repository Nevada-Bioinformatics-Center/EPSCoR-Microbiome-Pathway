/*
 * Module: KneadData
 * Description: This module performs kneading of metagenomic data.
 * Version: v1.29
 * Author: Kanishka Manna
 */

process RUN_MULTIQC {

    tag "MultiQC"

    label 'multiqc_conda'
    label 'multiqc_docker'

    publishDir "${params.output}/multiqc_out", mode: 'copy'

    input:
        path '*'

    output:
        path "multiqc_report.html", emit: multiqc_report
        path "multiqc_data"       , emit: multiqc_data

    script:
    """
    multiqc . --outdir ./ --force
    """
}