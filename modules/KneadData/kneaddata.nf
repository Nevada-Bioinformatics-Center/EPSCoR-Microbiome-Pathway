/*
 * Module: KneadData
 * Description: This module performs kneading of metagenomic data.
 * Version: v0.12.2
 * Author: Kanishka Manna and Hans Vasquez-Gross
 */

process KNEADING_DATA {
    tag "KneadData: ${sample_id}"

    label 'kneaddata'

    publishDir "${params.output}/kneaddata_out", mode: 'symlink'

    input:
        tuple val(sample_id), path(read), path(kneaddata_db_paths) 

    output:
        path( "${sample_id}.fastq.gz" ), emit: kneaddata_fastq
        path( "*.log" ), emit: kneaddata_log
        path( "fastqc/*_fastqc.html" ), emit: kneaddata_fastqc_html
        path( "fastqc/*_fastqc.zip" ), emit: kneaddata_fastqc_zip

    script:
        """
        kneaddata \\
            --input1 ${read[0]} \\
            --input2 ${read[1]} \\
            --threads ${params.kneaddata_threads} \\
            --processes ${params.kneaddata_processes} \\
            --reference-db ${kneaddata_db_paths} \\
            --output  "./" \\
            --output-prefix ${sample_id} \\
            --run-fastqc-start \\
            --run-fastqc-end \\
            --remove-intermediate-output \\
            --cat-final-output \\
            --run-trim-repetitive \\
            --sequencer-source ${params.kneaddata_sequencer_source} \\
            ${params.kneaddata_extra} \\
            --log ${sample_id}_kneaddata.log

        # gzip all FASTQ files
        gzip -f *.fastq
        """
}
