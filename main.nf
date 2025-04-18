#! /usr/bin/env nextflow

// To use DSL 2 will need to include this
nextflow.enable.dsl = 2

/*
 * ==================================================================================
 * EPSCoR-NASA Microbiome Pipeline: Pathway-level, consensus analysis of 
 * microbiome profiling data. This involves identifying and comparing 
 * metabolic pathways across microbial communities to derive a consensus 
 * understanding of their functional potential.
 * ==================================================================================
 */

// -------------------------------------------------------------
// `main.nf` is the pipeline script for this nextflow pipeline
// Should contain the following sections:
//      - Process definitions
//      - Pipeline structure
//      - Pipeline summary logs
// -------------------------------------------------------------


/*
 * Import modules
 * -------------------------------------------------------------
 */
// Import processes or subworkflows to be run in the pipeline
// Each of these is a separate .nf script saved in modules/ directory
include { QUALITY_CONTROL } from './modules/FastQC/fastqc.nf'
include { KNEADING_DATA } from './modules/KneadData/kneaddata.nf'
//include { TAXONOMIC_PROFILING } from './modules/MetaPhlAn/metaphlan.nf'
//include { FUNCTIONAL_PROFILING } from './modules/HUMAnN/humann.nf'


/*
 * Pipeline structure
 * -------------------------------------------------------------
 */
workflow  {

    /*
     * Input channels & parameters
     */

     // Show help message if the user specifies the `--help` flag
    if ( params.help ) {
        // Invoke the function above which prints the help message
        log.info params.help_message
        // Exit out and do not run anything else
        exit 1
    }

    // If the library input was not specified
    if ( !params.library == "PE" ) {
        exit 1, "Unsupported library preparation ${params.library}!"
    }

    // If the samplesheet input was specified
    if ( params.library == "PE") {

        samples_ch = Channel
                .fromPath( params.samplesheet )
                .ifEmpty { throw new IllegalArgumentException("Cannot find any sample-sheet! Please provide one by the input --samplesheet") }
                .splitCsv( header: ['sampleID', 'read1', 'read2'], sep: ',', skip: 1)
                .map{ row -> tuple([row.sampleID, [row.read1, row.read2]]) }
    }

    samples_ch.view() // check the paired ends to test the output, remove this once the testing is completed


    // Make a channel with the reference database specified by the user
    Channel.fromPath( params.kneaddata_db, checkIfExists: true ).set { kneaddata_DB }

    //Channel.fromPath( params.metaphlan_db, checkIfExists: true ).set { metaphlan_DB }

    //Channel.fromPath( params.humann_db, checkIfExists: true ).set { humann_DB }


    /*
     * Pipeline processes
     */
    QUALITY_CONTROL ( samples_ch )
    KNEADING_DATA ( samples_ch , kneaddata_DB)
    //TAXONOMIC_PROFILING ( KNEADING_DATA.out.fastq, metaphlan_DB )
    //FUNCTIONAL_PROFILING ( KNEADING_DATA.out.fastq, humann_DB )


    /*
     * Print pipeline execution summary
     */

}