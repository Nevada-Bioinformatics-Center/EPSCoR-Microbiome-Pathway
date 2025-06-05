
process RUN_MULTIQC {

    tag "MultiQC"

    label 'multiqc_conda'
    label 'multiqc_docker'
    label 'medium'

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