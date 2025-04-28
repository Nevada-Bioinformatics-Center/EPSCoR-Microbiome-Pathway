#! /usr/bin/env nextflow
//import java.nio.file.Paths

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

include { KNEADING_DATA } from './modules/KneadData/kneaddata.nf'
include { download_kneaddata_db } from './modules/KneadData/kneaddata_db.nf'
include { TAXONOMIC_PROFILING } from './modules/MetaPhlAn/metaphlan.nf'
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
        exit 0
    }

    // Check if samplesheet is specified
    if (!file(params.samplesheet).exists()) {
        exit 1, "Cannot find sample-sheet at ${params.samplesheet}! Please provide a --samplesheet!"
    }

    // If the library input was not specified
    if ( !params.library == "PE" ) {
        exit 1, "Unsupported library preparation ${params.library}!"
    }

    if (!params.kneaddata_db) {
        exit 1, "Please provide at least one --kneaddata_db in the format 'name:build'"
    }    


    // If the samplesheet input was specified
    if ( params.library == "PE") {

        samples_ch = Channel
            .fromPath("file:${params.samplesheet}", checkIfExists: true)
            .splitCsv(header: ['sampleID', 'read1', 'read2'], sep: ',', skip: 1)
            .map { row -> tuple(row.sampleID, [file(row.read1), file(row.read2)]) }

    }

    //samples_ch.view() // check the paired ends to test the output, remove this once the testing is completed


    // Make a channel with the reference database specified by the user

    //Channel.fromPath( params.metaphlan_db, checkIfExists: true ).set { metaphlan_DB }

    //Channel.fromPath( params.humann_db, checkIfExists: true ).set { humann_DB }

    def valid_db_options = [
        "human_genome"       : ["bowtie2", "bmtagger"],
        "human_transcriptome": ["bowtie2"],
        "ribosomal_RNA"      : ["bowtie2"],
        "mouse_C57BL"        : ["bowtie2"],
        "dog_genome"         : ["bowtie2"],
        "cat_genome"         : ["bowtie2"]
    ]

    def db_combinations = params.kneaddata_db
        .split(',')
        .collect { it.tokenize(':') } // Split each into [db, build]

    // Validate each combination
    db_combinations.each { entry ->
        if (entry.size() != 2) {
            exit 1, "❌ Invalid kneaddata_db entry: '${entry.join(':')}'. Must be in format: <db>:<build>"
        }
        def (db, build) = entry
        if (!valid_db_options.containsKey(db)) {
            exit 1, "❌ Unknown database: '${db}'. Must be one of: ${valid_db_options.keySet().join(', ')}"
        }
        if (!(build in valid_db_options[db])) {
            exit 1, "❌ Invalid build '${build}' for database '${db}'. Allowed builds: ${valid_db_options[db].join(', ')}"
        }
    }

    Channel
        .fromList(db_combinations)
        .map { db, build -> 
            def outdir = "${params.kneaddata_path}/${db}_${build}".replaceAll(/\/+/, '/')
            tuple(db, build, outdir)
        }
        .set { db_inputs }



    /*
     * Pipeline processes
     */

    // --------------------------------------------
    // Step 1: Quality control and processing FASTQ 
    //         paired-end reads with KneadData
    // ---------------------------------------------

    // Run download
    download_kneaddata_db(db_inputs)
        .collect()
        .map { db_list -> [ db_list ] }
        .set { kneaddata_db_ch }


    // Combine reads + database
    samples_ch
        .combine(kneaddata_db_ch)
        .set { kneading_inputs }

    //kneading_inputs.view()

    kneading_inputs
        .map { sample_id, read, dbs ->
            def out1 = java.nio.file.Paths.get("${params.output}/kneaddata_out/${sample_id}_paired_1.fastq").toFile()
            def out2 = java.nio.file.Paths.get("${params.output}/kneaddata_out/${sample_id}_paired_2.fastq").toFile()
            def exists = out1.exists() && out2.exists()

            if (exists && !params.force) {
                log.info "⏭️ Skipping ${sample_id} — output exists and --force not set"
            } else if (exists && params.force) {
                log.info "♻️ Rerunning ${sample_id} — output exists but --force is set"
            } else {
                log.info "✅ Will run ${sample_id} — output doesn't exist"
            }

            tuple(sample_id, read, dbs, exists)
        }
        .filter { sample_id, read, dbs, exists -> 
            return params.force || !exists
        }
        .map { sample_id, read, dbs, exists-> tuple(sample_id, read, dbs) }
        .set { kneading_inputs_filtered }

    kneading_inputs_filtered.view { it -> "🧪 KNEADDATA sample input: $it" }

    // Run KneadData
    KNEADING_DATA(kneading_inputs_filtered)

    KNEADING_DATA.out.kneaddata_fastq
                    .ifEmpty { exit 1, "❌ No output from KNEADING_DATA. Please check the input data or process configuration." }
                    .collectFile(name: { it.sample_id + '.fastq' })
                    .map { file ->
                        def sample_id = file.getName().replace('.fastq', '')
                        tuple (sample_id, file)
                     }
                    .set { merged_reads }


    // -----------------------------------------
    // Step 2: Profiling Taxonomy with MetaPhlAn
    // -----------------------------------------
    merged_reads.view{ "🧪 Merged reads for Taxonomy profiling: $it" }
    TAXONOMIC_PROFILING ( merged_reads )


    // ----------------------------------------
    // Step 3: Functional Profiling with HUMAnN
    // ----------------------------------------
    merged_reads.view{ "🧪 Merged reads for Functional profiling: $it" }
    //FUNCTIONAL_PROFILING ( KNEADING_DATA.out.fastq, humann_nucleotide_db, humann_protein_db, TAXONOMIC_PROFILING.out.profiled_taxa )


    /*
     * Pipeline event handler
     */
}



// End of the Pipeline. Goodbye!