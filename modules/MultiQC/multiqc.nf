/*
 * Module: MultiQC
 * Description: This module aggregates all QC and other info from Read Processing/Kneading Data into a report.
 * Version: v1.29
 * Author: Kanishka Manna and Hans Vasquez-Gross
 */

process RUN_MULTIQC_PROCESSED {

    tag "MultiQC"

    label 'multiqc'

    publishDir "${params.output}/multiqc_out/processed", mode: 'copy'

    input:
        path processed_fastqc

    output:
        path "multiqc_report.html", emit: multiqc_processed_report
        path "multiqc_data"

    script:
    """
    multiqc . --outdir ./ --force
    """
}


process RUN_MULTIQC_RAW {

    tag "MultiQC"

    label 'multiqc'

    publishDir "${params.output}/multiqc_out/raw", mode: 'copy'

    input:
        path raw_fastqc

    output:
        path "multiqc_report.html", emit: multiqc_raw_report
        path "multiqc_data"

    script:
    """
    multiqc . --outdir ./ --force
    """
}
