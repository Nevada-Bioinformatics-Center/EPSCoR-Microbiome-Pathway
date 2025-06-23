/*
 * Module: CPA
 * Description: This module performs consensus pathway analysis on metagenomic data.
 * Version: v0.1
 * Author: Kanishka Manna
 * Date: 2025-06-11
 */

process GENERATE_GOTERMS {
    tag "Gene Ontology Terms"

    label 'go_term_conda'
    label 'medium'

    output:
        path ("goterms/GOTerms.rds"), emit: goterms
    
    script:
    """
    mkdir -p goterms

    # Download files if not present
    [ -f goterms/gene2go.gz ] || wget -O goterms/gene2go.gz https://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz
    [ -f goterms/All_Data.gene_info.gz ] || wget -O goterms/All_Data.gene_info.gz https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/All_Data.gene_info.gz
    [ -f goterms/go.obo ] || wget -O goterms/go.obo http://purl.obolibrary.org/obo/go.obo

    if [ ! -f goterms/GOTerms.rds ]; then
        echo "Generating GOTerms.rds in goterms"
        goterms.R goterms/GOTerms.rds goterms/gene2go.gz goterms/All_Data.gene_info.gz goterms/go.obo
    else
        echo "GOTerms.rds already exists in goterms, skipping generation."
    fi
    """
}



process CPA_ANALYSIS {
    tag "Consensus Pathway Analysis"

    label 'cpa_conda'
    label 'high'
    
    publishDir "${params.outdir}/cpa_out", mode: 'copy'

    input:
        path (samplesheet)
        path (humann_dir)
        path (goterms)
 
    
    output:
        path ("/*")
    
    script:
    """
    cpa.R ${samplesheet} ${humann_dir} ${goterms}
    """
}