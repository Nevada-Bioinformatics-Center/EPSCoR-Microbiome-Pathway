/*
 * Module: CPA
 * Description: This module performs consensus pathway analysis on metagenomic data.
 * Version: v0.1
 * Author: Kanishka Manna
 */

// Generate GO terms for CPA analysis
process GENERATE_GOTERMS {
    tag "Gene Ontology Terms"

    label 'go_term'
    label 'medium'

    publishDir "${params.goterm_db}", mode: 'link', overwrite: true

    output:
        path ("GOTerms.rds"), emit: goterms
        path ("gene2go.gz")
        path ("All_Data.gene_info.gz")
        path ("go.obo")
        path (".done")
    
    when:
        ! file("${params.goterm_db}/GOTerms.rds").exists()
    
    script:
    """
    echo "Downloading and generating GO term resources..."
    
    # Download required files
    wget -O gene2go.gz https://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz
    wget -O All_Data.gene_info.gz https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/All_Data.gene_info.gz
    wget -O go.obo http://purl.obolibrary.org/obo/go.obo

    # Run the script using the dynamic path
    ${task.ext.script_path} GOTerms.rds gene2go.gz All_Data.gene_info.gz go.obo

    # Create the done file to signal completion
    touch .done
    """
}



// Process to extract metadata from the samplesheet
process EXTRACT_METAINFO {
    tag "Extract Meta information"

    label 'cpa'
    label 'low'

    input:
        path (samplesheet)
    
    output:
        path ("cpa_metainfo.tsv"), emit: metainfo
    
    script:
    """
    cut -d, -f 1,2 ${samplesheet} > cpa_metainfo.tsv
    """
}



// Process to run CPA analysis
process CONSENSUS_PATHWAY_ANALYSIS {
    tag "Consensus Pathway Analysis"

    label 'cpa'
    label 'high'
    
    publishDir "${params.output}", mode: 'move'

    input:
        path (samplesheet)
        path (humann_dir)
        path (goterms)
 
    
    output:
        path ("*")
    
    script:
    """
    ${task.ext.script_path} ${samplesheet} ${humann_dir} ${goterms}
    """
}
