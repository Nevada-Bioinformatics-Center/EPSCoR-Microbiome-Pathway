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

include { KNEADING_DATA } from './modules/KneadData/kneaddata.nf'
include { TAXONOMIC_PROFILING  } from './modules/MetaPhlAn/metaphlan.nf'
include { MERGE_TAXONOMIC_TABLES } from './modules/MetaPhlAn/metaphlan_utility.nf'
include { TAXONOMIC_VISUALIZATION } from './modules/MetaPhlAn/metaphlan_utility.nf'
include { FUNCTIONAL_PROFILING } from './modules/HUMAnN/humann.nf'
include { NORMALIZE_PATHWAY_ABUNDANCE } from './modules/HUMAnN/humann_utility.nf'
include { JOIN_PATHWAY_ABUNDANCE } from './modules/HUMAnN/humann_utility.nf'
include { DESCRIPTIVE_PROFILING } from './modules/HUMAnN/humann_utility.nf'
include { GENERATE_GOTERMS } from './modules/CPA/cpa.nf'
include { EXTRACT_METAINFO } from './modules/CPA/cpa.nf'
include { CONSENSUS_PATHWAY_ANALYSIS } from './modules/CPA/cpa.nf'
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
            .fromPath("${params.samplesheet}", checkIfExists: true)
            .splitCsv(header: true, sep: ',', skip: 0)
            .map { row -> 
            if (!row.fastq_1 || !row.fastq_2) {
                exit 1, "Each row must contain 'sample', 'fastq_1', and 'fastq_2'." 
            }
            tuple(row.sample, row.factors, [file(row.fastq_1), file(row.fastq_2)])
            }
    
    // View parsed samples
    samples_ch.view { it -> "Parsed sample: $it" }
        
    // Check if the samplesheet has 'exp_conditions'
    def samplesheet_header = file(params.samplesheet).text.readLines().first().split(',')*.trim()
    def exp_conditions_idx = samplesheet_header.indexOf('exp_conditions')
    def has_exp_conditions = exp_conditions_idx != -1
//        def samplesheet_header = file(params.samplesheet).text.readLines().first().split(',')
//        def has_exp_conditions = samplesheet_header*.trim().contains('exp_conditions')

        // Get all experiment conditions from the samplesheet
//        def exp_conditions_values = file(params.samplesheet)
//            .text
//            .readLines()
//            .drop(1) // skip header
//            .collect { it.split(',')[samplesheet_header.indexOf('exp_conditions')].trim() }
//            .unique()

        def exp_conditions_values = has_exp_conditions
            ? file(params.samplesheet)
                .text
                .readLines()
                .drop(1) // skip header
                .collect { it.split(',')[exp_conditions_idx].trim() }
                .unique()
            : []
        
        // Check if exp_conditions has at least two unique values
        def exp_conditions_ok = has_exp_conditions && exp_conditions_values.size() > 1


    /*
        -----------------------------------------------------
        Step 1: Quality Control & Processing paired-end FASTQ
        -----------------------------------------------------
    */

    // Check if kneaddata_db parameter is provided
    if (! params.kneaddata_db) {
        exit 1, "Please provide Kneaddata databases using the --kneaddata_db"
    }

    // Channel for KneadData database
    Channel.fromPath(params.kneaddata_db, checkIfExists: true)
        .set { kneaddata_db_ch }

    // Combine samples with the database for KneadData
    samples_ch
        .combine(kneaddata_db_ch)
        .map { sample_id, _factors, read, dbs -> tuple(sample_id, read, dbs) }
        .set { kneaddata_inputs }

    
    // Run KneadData
    KNEADING_DATA(kneaddata_inputs)



    /*
        --------------------------
        Step 2: Taxonomy Profiling
        --------------------------
    */

    // Check if MetaPhlAn database is provided
    if (!params.metaphlan_db) {
        exit 1, "Please provide MetaPhlAn database using the --metaphlan_db"
    }

    // Channel for MetaPhlAn database
    Channel.fromPath(params.metaphlan_db, checkIfExists: true)
        .set { metaphlan_db_dir_ch }

    // Combine merged reads with MetaPhlAn database path
    KNEADING_DATA.out.kneaddata_fastq
        .collectFile()
        .map { fastq_file -> 
            def sample_id = fastq_file.getBaseName().replace('.fastq.gz', '').replace('.fastq', '')
            tuple(sample_id, fastq_file)
        }
        .combine(metaphlan_db_dir_ch)
        .map { sample_id, reads, db_dir -> 
            tuple(sample_id, reads, db_dir) 
        }
       .set { metaphlan_inputs }

    metaphlan_inputs.view { "MetaPhlAn inputs: $it" }

    // Run MetaPhlAn for taxonomic profiling
    TAXONOMIC_PROFILING(metaphlan_inputs)


    // Collect all MetaPhlAn profile.tsv files
    TAXONOMIC_PROFILING.out.profiled_taxa
        .collect()
        .set { metaphlan_tables_ch }

    Channel.fromPath(params.samplesheet, checkIfExists: true)
        .set { samplesheet_file_ch }
    

    if (params.metaphlan_extra_analysis) {
        // Run MetaPhlAn utility to merge tables
        MERGE_TAXONOMIC_TABLES(metaphlan_tables_ch)

        MERGE_TAXONOMIC_TABLES.out.merged_taxa_profile.set { merged_taxa_profile_ch }

        // Run MetaPhlAn visualization
        TAXONOMIC_VISUALIZATION(
            merged_taxa_profile_ch,
            samplesheet_file_ch
        )
    }



    /*
        ----------------------------
        Step 3: Functional Profiling
        ----------------------------
    */

    // Check if HUMAnN nucleotide database is provided
    if (!params.humann_nucleotide_db) {
        exit 1, "Please provide HUMAnN nucleotide database using --humann_nucleotide_db"
    }

    // Channel for HUMAnN nucleotide database
    Channel.fromPath(params.humann_nucleotide_db)
        .set { humann_nucleotide_db_ch }
    
    // Check if HUMAnN protein database is provided
    if (!params.humann_protein_db) {
        exit 1, "Please provide HUMAnN protein database using --humann_protein_db"
    }

    // Channel for HUMAnN protein database
    Channel.fromPath(params.humann_protein_db)
        .set { humann_protein_db_ch }

    // Channel for profiled taxa from MetaPhlAn
    TAXONOMIC_PROFILING.out.profiled_taxa
        .collectFile()
        .map { taxa_file ->
            def sample_id = taxa_file.getName().replaceFirst(/_profile\.tsv$/,'')
            tuple(sample_id, taxa_file)
        }
        .set { sample_taxa_ch }

    // Combine files for Functional profiling inputs
    KNEADING_DATA.out.kneaddata_fastq
        .collectFile()
        .map { fastq_file -> 
            def sample_id = fastq_file.getBaseName().replace('.fastq.gz', '').replace('.fastq', '')
            tuple(sample_id, fastq_file)
        }
        .set { sample_fastq_ch }

    // Combine all inputs for HUMAnN
    sample_fastq_ch
    .combine(humann_nucleotide_db_ch)
        .map { sample_id, reads, nuc_db -> tuple(sample_id, reads, nuc_db) }
    .combine(humann_protein_db_ch)
        .map { sample_id, reads, nuc_db, prot_db -> tuple(sample_id, reads, nuc_db, prot_db) }
    .join(sample_taxa_ch)
        .map { sample_id, reads, nuc_db, prot_db, taxa ->
            tuple(sample_id, reads, nuc_db, prot_db, taxa)
        }
        .set { humann_inputs }

    humann_inputs.view { "HUMAnN inputs: $it" }

    // Run HUMAnN functional profiling
    FUNCTIONAL_PROFILING(humann_inputs)



   if (exp_conditions_ok) {

    /*
        ----------------------------
        Step 4: CPA analysis
        ----------------------------
    */

    // Generate GO terms from pathway data
    GENERATE_GOTERMS()

    // Channel for existing GO terms
    ext_goterms = Channel.fromPath("${params.goterm_db}/GOTerms.rds")

    // Merge the generated GO terms with existing ones
    goterms_ch = GENERATE_GOTERMS.out.goterms.mix(ext_goterms).first()



    // Collect HUMAnN gene family output files into a directory for CPA input
    FUNCTIONAL_PROFILING.out.gene_fam
    .collect()
    .map { files ->
        def outdir = file("humann_genefam_dir")
        outdir.mkdirs()
        files.each { f -> 
            def dst = outdir.resolve(f.getName())
            f.copyTo(dst)
        }
        return outdir
    }
    .set { humann_genefam_ch }
    
    humann_genefam_ch.view { "Directory prepared for CPA: $it" }

    // Extract sample metadata information from the samplesheet
    Channel.fromPath(params.samplesheet, checkIfExists: true)
            .set { samplesheet_file_ch }

    EXTRACT_METAINFO(samplesheet_file_ch)

    // Run CPA analysis
   CONSENSUS_PATHWAY_ANALYSIS(EXTRACT_METAINFO.out.metainfo, humann_genefam_ch, goterms_ch)
   }
   


    /*
        -------------------------------------
        Step 5: Descriptive Pathway Profiling
        -------------------------------------
    */

    // Collect HUMAnN pathway abundance output files and pair with sample IDs
    FUNCTIONAL_PROFILING.out.path_abundance
        .collectFile()
        .map { pathabundance_file -> 
            def sample_id = pathabundance_file.getBaseName().replace('_pathabundance.tsv', '')
            tuple(sample_id, pathabundance_file)
        }
        .set { pathabundance_ch }

    // Normalize pathway abundance tables
    NORMALIZE_PATHWAY_ABUNDANCE (pathabundance_ch)

    // Collect normalized pathway abundance files into a directory
    NORMALIZE_PATHWAY_ABUNDANCE.out.renorm_pathabundance
        .collect()
        .map { files ->
            def outdir = file("humann_renorm_dir")
            outdir.mkdirs()
            files.each { f -> 
                def dst = outdir.resolve(f.getName())
                f.copyTo(dst)
            }
            return outdir
        }
        .set { renorm_pathabun_ch }

    // Join normalized pathway abundance tables
    JOIN_PATHWAY_ABUNDANCE(renorm_pathabun_ch)

    // Run descriptive profiling
    DESCRIPTIVE_PROFILING(JOIN_PATHWAY_ABUNDANCE.out.join_table)



    /*
        ----------------------------
        Step 6: MultiQC reports
        ----------------------------
    */

    // Run multiqc
    RUN_MULTIQC(
        KNEADING_DATA.out.kneaddata_fastqc_html.mix(
            KNEADING_DATA.out.kneaddata_fastqc_zip,
            TAXONOMIC_PROFILING.out.profiled_taxa_txt
        ).collect()
    )
}


// End of the Pipeline. Goodbye!