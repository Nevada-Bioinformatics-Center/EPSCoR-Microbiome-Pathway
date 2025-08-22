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

    publishDir "${params.goterm_db}", mode: 'symlink'

    output:
        path ("GOTerms.rds"), emit: goterms
        path ("gene2go.gz")
        path ("All_Data.gene_info.gz")
        path ("go.obo")
        path (".done")
    
    when:
        ! file("${params.goterm_db}/GOTerms.rds").exists() || 
        ! file("${params.goterm_db}/gene2go.gz").exists() || 
        ! file("${params.goterm_db}/All_Data.gene_info.gz").exists() || 
        ! file("${params.goterm_db}/go.obo").exists() ||
        ! file("${params.goterm_db}/.done").exists()
    
    script:
    """
    db_dir="."
    
    if [ ! -f "\$db_dir/.done" ]; then
        echo "Downloading and generating GO term resources in \$db_dir"
        [ -f gene2go.gz ] || wget -O gene2go.gz https://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz
        [ -f All_Data.gene_info.gz ] || wget -O All_Data.gene_info.gz https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/All_Data.gene_info.gz
        [ -f go.obo ] || wget -O go.obo http://purl.obolibrary.org/obo/go.obo

        goterms.R GOTerms.rds gene2go.gz All_Data.gene_info.gz go.obo
        touch .done
    else
        echo "GOTerms.rds already exists in \$db_dir, skipping download and generation."
    fi
    """
}



// Process to extract metadata from the samplesheet
process EXTRACT_METAINFO {
    tag "Extract Meta information"

    label 'cpa_conda'
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

    label 'cpa_conda'
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
    cpa.R ${samplesheet} ${humann_dir} ${goterms}
    """
}
