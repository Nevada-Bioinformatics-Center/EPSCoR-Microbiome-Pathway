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
include { TAXONOMIC_PROFILING  } from './modules/MetaPhlAn/metaphlan.nf'
include { FUNCTIONAL_PROFILING } from './modules/HUMAnN/humann.nf'
include { DOWNLOAD_METPHLAN_DB } from './modules/MetaPhlAn/metaphlan.nf'
//include { COPY_METAPHLAN_DB } from './modules/MetaPhlAn/metaphlan.nf'




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
    Channel.fromPath( params.nucleotide_db, checkIfExists: true ).set { nucleotide_db_ch }
    Channel.fromPath( params.protein_db, checkIfExists: true ).set { protein_db_ch }

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
        .map { sample_id, read, dbs ->
            def out1 = file("${params.output}/kneaddata_out/${sample_id}_paired_1.fastq")
            def out2 = file("${params.output}/kneaddata_out/${sample_id}_paired_2.fastq")

            def exists = out1.exists() && out2.exists()
            tuple(sample_id, read, dbs, exists)
        }
        .set { kneaddata_status_ch }

    //kneading_inputs.view()


    kneaddata_status_ch
        .filter { sample_id, read, dbs, exists -> !exists || params.force }
        .map    { sample_id, read, dbs, exists -> tuple(sample_id, read, dbs) }
        .set    { kneaddata_inputs_filtered }

    kneaddata_inputs_filtered.view { it -> "🧪 KNEADDATA sample input: $it" }

    // Run KneadData
    KNEADING_DATA(kneaddata_inputs_filtered)


    kneaddata_status_ch
        .filter { sample_id, read, dbs, exists -> exists && !params.force }
        .map { sample_id, read, dbs, exists ->
            def out1 = file("${params.output}/kneaddata_out/${sample_id}_paired_1.fastq")
            def out2 = file("${params.output}/kneaddata_out/${sample_id}_paired_2.fastq")
            tuple(sample_id, [out1, out2])
        }
        .concat(KNEADING_DATA.out.kneaddata_fastq)
        .set { merged_reads }



    // -----------------------------------------
    // Step 2: Profiling Taxonomy with MetaPhlAn
    // -----------------------------------------
    // Create channel for MetaPhlAn DB path
    Channel.value(params.metaphlan_dbpath).set { metaphlan_db_dir_ch }

    // Download MetaPhlAn DB
    DOWNLOAD_METPHLAN_DB(metaphlan_db_dir_ch)

    DOWNLOAD_METPHLAN_DB.out.metaphlan_db_dir
        .set { metaphlan_db_dir_out_ch }


    // Combine sample reads with the MetaPhlAn database path
    //merged_reads
    //    .combine(DOWNLOAD_METPHLAN_DB.out.metaphlan_db)
    //    .map { sample_id, reads, dbpath -> tuple(sample_id, reads, dbpath) }
    //    .set { metaphlan_inputs }
    merged_reads
    .combine(metaphlan_db_dir_out_ch)
    .map { sample_id, reads, db_dir -> tuple(sample_id, reads, db_dir) }
    .set { metaphlan_inputs }

    metaphlan_inputs.view { "🧪 MetaPhlAn input: $it" }
    TAXONOMIC_PROFILING(metaphlan_inputs)

    //TAXONOMIC_PROFILING.out.profiled_taxa
    //.map { file -> 
    //    def sample_id = file.getName().replace('_profile.tsv', '') // Extract sample_id from file name
    //    tuple(sample_id, file)
    //}
    //.set { profiled_taxa_mapped }


    // ----------------------------------------
    // Step 3: Functional Profiling with HUMAnN
    // ----------------------------------------
    
    merged_reads.combine(TAXONOMIC_PROFILING.out.profiled_taxa)
                .map { sample_id, reads, profiled_taxa -> 
                    tuple(sample_id, reads, profiled_taxa)
                }
    .set { functional_inputs }
    //merged_reads.combine(profiled_taxa_mapped)
    //            .filter { merged, profiled_taxa -> 
    //                merged[0] == profiled_taxa[0] // Match sample_id
    //            }
    //            .map { merged, profiled_taxa -> 
    //                tuple(merged[0], merged[1], profiled_taxa[1]) // sample_id, reads, taxonomic_profile
    //            }
    //            .set { functional_inputs }
    functional_inputs.view { "🧪 Functional profiling input: $it" }

    FUNCTIONAL_PROFILING ( functional_inputs, nucleotide_db_ch, protein_db_ch )


    /*
     * Pipeline event handler
     */
}



// End of the Pipeline. Goodbye!