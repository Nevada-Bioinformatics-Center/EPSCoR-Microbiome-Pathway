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
// Contains the following sections:
//      - Import modules
//      - Pipeline structure
//      - Pipeline summary logs
// -------------------------------------------------------------


/*
 * Import modules
 * -------------------------------------------------------------
 */
// Import processes or subworkflows to be run in the pipeline
// Each of these is a separate .nf script saved in modules/ directory

include { DOWNLOAD_KNEADDATA_DB } from './modules/KneadData/kneaddata_db.nf'
include { KNEADING_DATA } from './modules/KneadData/kneaddata.nf'
include { DOWNLOAD_METAPHLAN_DB } from './modules/MetaPhlAn/metaphlan_db.nf'
include { TAXONOMIC_PROFILING  } from './modules/MetaPhlAn/metaphlan.nf'
include { DOWNLOAD_HUMANN_NUCLEOTIDE_DB } from './modules/HUMAnN/humann_db.nf'
include { DOWNLOAD_HUMANN_PROTEIN_DB } from './modules/HUMAnN/humann_db.nf'
include { FUNCTIONAL_PROFILING } from './modules/HUMAnN/humann.nf'
include { NORMALIZE_PATHWAY_ABUNDANCE } from './modules/HUMAnN/humann_utility.nf'
include { JOIN_PATHWAY_ABUNDANCE } from './modules/HUMAnN/humann_utility.nf'
include { DESC_PROFILING } from './modules/HUMAnN/humann_utility.nf'
include { GENERATE_GOTERMS } from './modules/CPA/cpa.nf'
include { EXTRACT_METAINFO } from './modules/CPA/cpa.nf'
include { CPA_ANALYSIS } from './modules/CPA/cpa.nf'
include { RUN_MULTIQC } from './modules/MultiQC/multiqc.nf'





/*
 * Pipeline structure
 * -------------------------------------------------------------
 */
workflow  {

    /*
     *  Pipeline Logic
     */
    
    /*
        Show help message if the user specifies the `--help` flag
    */
    if ( params.help ) {
        // Invoke the function above which prints the help message
        log.info params.help_message
        // Exit out and do not run anything else
        exit 0
    }

    /*
        ---------------------------------
        Samplesheet: Parsing & Validation
        ---------------------------------
    */
    // Check if samplesheet is specified
    if (!file(params.samplesheet).exists()) {
        exit 1, "Cannot find sample-sheet at ${params.samplesheet}! Please provide a --samplesheet!"
    }    

    // Parse the samplesheet and validate paired-end format
        samples_ch = Channel
            .fromPath("file:${params.samplesheet}", checkIfExists: true)
            .splitCsv(header: true, sep: ',', skip: 0)
            .map { row -> 
            if (!row.fastq_1 || !row.fastq_2) {
                exit 1, "Each row must contain 'sample', 'factors', 'fastq_1', and 'fastq_2'." 
            }
            tuple(row.sample, row.factors, [file(row.fastq_1), file(row.fastq_2)])
            }
        
//        def samplesheet_header = file(params.samplesheet).text.readLines().first().split(',')
//        def has_factor = samplesheet_header*.trim().contains('factor')


    /*
        -----------------------------------------------------
        Step 1: Quality Control & Processing paired-end FASTQ
        -----------------------------------------------------
    */

    // Check if `--kneaddata_db` parameter is provided
    if (!params.kneaddata_db) {
        exit 1, "Please provide at least one --kneaddata_db in the format 'name:build'"
    }

    // Valid Kneaddata database options
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

    // Validate each combinations
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

// Define DB combinations with expected paths
Channel
    .fromList(db_combinations)
    .map { db, build ->
        def db_dir = "${params.database}/kneaddata_db/${db}_${build}"
        def done_file = file("${db_dir}")
        tuple(db, build, db_dir, done_file.exists())
    }
    .set { db_inputs_with_check }

// split into two channels: exists or not exists
db_exists_ch = db_inputs_with_check
    .filter { _db, _build, _db_dir, exists -> exists }
    .map { _db, _build, db_dir, _exists -> file(db_dir) }

db_not_exists_ch = db_inputs_with_check
    .filter { _db, _build, _db_dir, exists -> !exists }
    .map { db, build, _db_dir, _exists -> tuple(db, build) }
    | DOWNLOAD_KNEADDATA_DB

// Merge both channels (existing + downloaded)
kneaddata_db_ch = db_exists_ch.mix(db_not_exists_ch)
    .collect()
    .map { db_list -> [db_list] }


    // Combine reads + database
    samples_ch
        .combine(kneaddata_db_ch)
        .map { sample_id, factors, reads, dbs ->
            def out1 = file("${params.output}/kneaddata_out/${sample_id}_paired_1.fastq")
            def out2 = file("${params.output}/kneaddata_out/${sample_id}_paired_2.fastq")

            def exists = out1.exists() && out2.exists()
            tuple(sample_id, factors, reads, dbs, exists)
        }
        .set { kneaddata_status_ch }

    // Final check for Kneaddata inputs
    kneaddata_status_ch
        .filter { _sample_id, _factors, _reads, _dbs, exists -> !exists || params.force }
        .map    { sample_id, _factors, read, dbs, _exists -> tuple(sample_id, read, dbs) }
        .set    { kneaddata_inputs_filtered }

    kneaddata_inputs_filtered.view { it -> "🧪 KNEADDATA sample input: $it" }

    // Run KneadData
    KNEADING_DATA(kneaddata_inputs_filtered)
    


    /*
        --------------------------
        Step 2: Taxonomy Profiling
        --------------------------
    */

    // Download MetaPhlAn database
    DOWNLOAD_METAPHLAN_DB()

    // Emit the MetaPhlAn database directory as queue channel
    Channel.value("${params.database}/metaphlan_db")
        .set { metaphlan_db_dir_out_ch }

    metaphlan_db_dir_out_ch.view {"$it"}

    // Combine merged reads with MetaPhlAn database path
    KNEADING_DATA.out.kneaddata_fastq
        .collectFile()
        .map { fastq_file -> 
            def sample_id = fastq_file.getBaseName().replace('.fastq', '')
            tuple(sample_id, fastq_file)
        }
        .combine(metaphlan_db_dir_out_ch)
        .map { sample_id, reads, db_dir -> 
            tuple(sample_id, reads, db_dir) 
        }
        .set { metaphlan_inputs }

    metaphlan_inputs.view { "🧪 MetaPhlAn input: $it" }

    // Run MetaPhlAn for taxonomic profiling
    TAXONOMIC_PROFILING(metaphlan_inputs)

//    TAXONOMIC_PROFILING.out.profiled_taxa.set { profiled_taxa }



    /*
        ----------------------------
        Step 3: Functional Profiling
        ----------------------------
    */
    
    // Define valid database types for HUMAnN
//    def valid_humann_db_options = [
//        "chocophlan" : ["full"],
//        "uniref" : ["uniref50_diamond", "uniref90_diamond", "uniref50_ec_filtered_diamond", "uniref90_ec_filtered_diamond"]
//    ]

    // Validate and process the humann nucleotide database parameter
//    if (params.humann_nucleotide_db) {
//        def nucleotide_db_combo = params.humann_nucleotide_db
//            .split(',')
//            .collect { it.tokenize(':') } //split each entry [db, build]

//            nucleotide_db_combo.each { entry ->
//                if (entry.size() != 2) {
//                    exit 1, "❌ Invalid humann nucleotide database entry: '${entry.join(':')}'. Must be in format: <nuc_db>:<nuc_build>"
//                }
//                def (nuc_db, nuc_build) = entry
//                if (!valid_humann_db_options.containsKey(nuc_db)) {
//                    exit 1, "❌ Unknown database: '${nuc_db}'. Must be one of: ${valid_humann_db_options.keySet().join(', ')}"
//                }
//                if (!(nuc_build in valid_humann_db_options[nuc_db])) {
//                    exit 1, "❌ Invalid build '${nuc_build}' for database '${nuc_db}'. Allowed builds: ${valid_humann_db_options[nuc_db].join(', ')}"
//                }
//            }

//            Channel.fromList(nucleotide_db_combo)
//                .map { db, build -> tuple( db, build ) }
//                .set { nucleotide_db_inputs }

            // Download HUMAnN nucleotide database
//            DOWNLOAD_HUMANN_NUCLEOTIDE_DB(nucleotide_db_inputs)
//                .collect()
//                .map { db_list -> [db_list] }
//                .set { humann_nucleotide_db_ch }
//    }
//    else {
//        println "⚠️ HUMAnN nucleotide database not provided. Functional profiling may not work as expected."
//    }

//    humann_nucleotide_db_ch.view { "$it" }

    
    // Validate and process the humann protein database parameter
//    if (params.humann_protein_db) {
//        def protein_db_combo = params.humann_protein_db
//            .split(',')
//            .collect { it.tokenize(':') } //split each entry [db, build]

//            protein_db_combo.each { entry ->
//                if (entry.size() != 2) {
//                    exit 1, "❌ Invalid humann protein database entry: '${entry.join(':')}'. Must be in format: <prot_db>:<prot_build>"
//                }
//                def (prot_db, prot_build) = entry
//                if (!valid_humann_db_options.containsKey(prot_db)) {
//                    exit 1, "❌ Unknown database: '${prot_db}'. Must be one of: ${valid_humann_db_options.keySet().join(', ')}"
//                }
//                if (!(prot_build in valid_humann_db_options[prot_db])) {
//                    exit 1, "❌ Invalid build '${prot_build}' for database '${prot_db}'. Allowed builds: ${valid_humann_db_options[prot_db].join(', ')}"
//                }
//            }

//            Channel.fromList(protein_db_combo)
//                .map { db, build -> tuple( db, build ) }
//                .set { protein_db_inputs }

            // Download HUMAnN protien database
//            DOWNLOAD_HUMANN_PROTEIN_DB(protein_db_inputs)
//                .collect()
//                .map { db_list -> [db_list] }
//                .set { humann_protein_db_ch }
//    }
//    else {
//        println "⚠️ HUMAnN protein database not provided. Functional profiling may not work as expected."
//    }

//    humann_protein_db_ch.view { "$it" }


    // Collect KneadData output FASTQ files and pair with sample IDs
//    KNEADING_DATA.out.kneaddata_fastq
//        .collectFile()
//        .map { fastq_file -> 
//            def sample_id = fastq_file.getBaseName().replace('.fastq', '')
//            tuple(sample_id, fastq_file)
//        }
//        .set { sample_fastq_ch }

    // Collect profiled taxa files and pair with sample IDs
//    profiled_taxa
//        .collectFile()
//        .map { taxa_file ->
//            def sample_id = taxa_file.getName().replace('_profile.tsv', '')
//            tuple(sample_id, taxa_file)
//         }
//        .set { sample_taxa_ch }

    // Combine sample FASTQ with HUMAnN nucleotide DB
//    sample_fastq_ch
//        .combine(humann_nucleotide_db_ch)
//        .map { sample_id, reads, nuc_db ->
//            tuple(sample_id, reads, nuc_db)
//        }
        // Combine with HUMAnN protein DB
//        .combine(humann_protein_db_ch)
//        .map { sample_id, reads, nuc_db, prot_db ->
//            tuple(sample_id, reads, nuc_db, prot_db)
//        }
        // Join with profiled taxa for each sample
//        .join(sample_taxa_ch)
//        .map { sample_id, reads, nuc_db, prot_db, taxa ->
//            tuple(sample_id, reads, nuc_db, prot_db, taxa)
//        }
//        .set { humann_inputs }

    // View HUMAnN input tuples for debugging
//    humann_inputs.view()

    // Run HUMAnN functional profiling
//    FUNCTIONAL_PROFILING (humann_inputs)



//   if (has_factor) {

    /*
        ----------------------------
        Step 4: CPA analysis
        ----------------------------
    */

    // Generate GO terms from pathway data
//    GENERATE_GOTERMS()

    // Collect HUMAnN gene family output files into a directory for CPA input
//    FUNCTIONAL_PROFILING.out.gene_fam
//    .collect()
//    .map { files ->
//        def outdir = file("humann_genefam_dir")
//        outdir.mkdirs()
//        files.each { f -> 
//            def dst = outdir.resolve(f.getName())
//            f.copyTo(dst)
//        }
//        return outdir
//    }
//    .set { humann_genefam_ch }
//    humann_genefam_ch.view { "🧾 Directory prepared for CPA: $it" }

    // Extract sample metadata information from the samplesheet
//    Channel.fromPath(params.samplesheet, checkIfExists: true)
//            .set { samplesheet_file_ch }
//    EXTRACT_METAINFO(samplesheet_file_ch)

    // Run CPA analysis
//    CPA_ANALYSIS(EXTRACT_METAINFO.out.metainfo, humann_genefam_ch, GENERATE_GOTERMS.out.goterms)
//   }
//   else {

    /*
        -------------------------------------
        Step 5: Descriptive Pathway Profiling
        -------------------------------------
    */

    // Collect HUMAnN pathway abundance output files and pair with sample IDs
//    FUNCTIONAL_PROFILING.out.path_abundance
//        .collectFile()
//        .map { pathabundance_file -> 
//            def sample_id = pathabundance_file.getBaseName().replace('_pathabundance.tsv', '')
//            tuple(sample_id, pathabundance_file)
//        }
//        .set { pathabundance_ch }

    // Normalize pathway abundance tables
//    NORMALIZE_PATHWAY_ABUNDANCE (pathabundance_ch)

    // Collect normalized pathway abundance files into a directory
//    NORMALIZE_PATHWAY_ABUNDANCE.out.renorm_pathabundance
//        .collect()
//        .map { files ->
//            def outdir = file("humann_renorm_dir")
//            outdir.mkdirs()
//            files.each { f -> 
//                def dst = outdir.resolve(f.getName())
//                f.copyTo(dst)
//            }
//            return outdir
//        }
//        .set { renorm_pathabun_ch }

//    renorm_pathabun_ch.view()

    // Join normalized pathway abundance tables
//    JOIN_PATHWAY_ABUNDANCE(renorm_pathabun_ch)

    // Run descriptive profiling
//    DESC_PROFILING(JOIN_PATHWAY_ABUNDANCE.out.join_table)

//   }



    /*
        ----------------------------
        Step 6: MultiQC reports
        ----------------------------
    */

    // Run multiqc
//    RUN_MULTIQC(
//        KNEADING_DATA.out.kneaddata_fastqc_zip.mix(
//            KNEADING_DATA.out.kneaddata_fastq,
//            KNEADING_DATA.out.kneaddata_log,
//            KNEADING_DATA.out.kneaddata_fastqc_html,
//            TAXONOMIC_PROFILING.out.profiled_taxa,
//            TAXONOMIC_PROFILING.out.metaphlan_log,
//            FUNCTIONAL_PROFILING.out.gene_fam,
//            FUNCTIONAL_PROFILING.out.path_abundance,
//            FUNCTIONAL_PROFILING.out.path_coverage,
//            FUNCTIONAL_PROFILING.out.humann_log
//        ).collect()
//    )


    /*
     * Pipeline event handler
     */
}


// End of the Pipeline. Goodbye!