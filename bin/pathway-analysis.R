#library(tidyverse) # install.packages("tidyverse")
#library(RCPA) # devtools::install_github("tinnlab/RCPA") or install.packages("RCPA")
#library(fgsea) # BiocManager::install("fgsea")
#library(pbmcapply) # BiocManager::install("pbmcapply")
#library(stringr)
#library(UniProt.ws)

setwd("~/Documents/sandbox/projects/EPSCoR-NASA/Melanie_Hess_codes/")
exportDir = "export_renamed"
# dataset = "OSD-212_renamed"
dataset = "OSD-69_renamed"
dir.create(exportDir, showWarnings = FALSE)

library(pacman)
p_load(tidyverse, RCPA, fgsea,pbmcapply,stringr,UniProt.ws)

#source("/home/melaniehess/Documents/NASA_Microbiomes/utils.R")
source("/home/kmanna/Documents/sandbox/projects/EPSCoR-NASA/Melanie_Hess_codes/utils.R")

if (file.exists("GOTerms.rds")) {
  GOTerms <- readRDS("GOTerms.rds")
} else {
  GOTerms <- getGOTerms()
  saveRDS(GOTerms, "GOTerms.rds")
}

## Can be shortened A LOT
# load GLDS-212 abundance files
gaFiles <- list.files(paste("humann_out/OSD-69_raw_sub/",dataset,sep=""), pattern = "*_genefamilies.tsv", full.names = TRUE, recursive = TRUE)
abundances <- lapply(gaFiles, function(file) {
  df <- read_tsv(file,comment="#",col_names=F,show_col_types = FALSE)
  colnames(df) <- c("Uniprot",gsub("genefamilies.tsv","Abundance-RPKs",basename(file)))
  df
})
allGenes <- Reduce(union, lapply(abundances, `[[`, "Uniprot"))
abundanceMat <- lapply(abundances, function(df) {
  data.frame(
    Uniprot = allGenes
  ) %>%
    left_join(df, by = "Uniprot") %>%
    mutate(across(-Uniprot, ~replace_na(., 0))) %>%
    dplyr::select(-Uniprot)
}) %>% bind_cols()
rownames(abundanceMat) <- allGenes
abundanceMat <- abundanceMat[!str_detect(rownames(abundanceMat), "unclassified|UNMAPPED"),]
rownames(abundanceMat) <- str_replace(rownames(abundanceMat), "UniRef90_", "")

write.csv(abundanceMat,paste(exportDir,"/abundances.csv",sep=""),row.names=T,quote=F)

# map to gene symbol
mapping <- UniProt.ws::mapUniProt(from = "UniProtKB_AC-ID",
                                  to = "Gene_Name",
                                  columns = character(0L),
                                  query = rownames(abundanceMat),
                                  verbose = FALSE,
                                  debug = FALSE,
                                  paginate = FALSE)

## Some of them don't map:
#  UPI00010EE582, UPI0009B56B74, A0A287LSF3, U5EL12, UPI000D58D589, UPI0009337703, A6N2T4, P24650, UPI000DED975A, 
#  A0A0A9RWP6, UPI0000113354, B3WYW0, UPI000988B2A4, UPI0009AD3D0E, UPI0005CCFA93, A0A390Y369, UPI000C7127F2, 
#  A0A390TGH2, UPI0003FBC9D6, A0A387JN55, A4QMB1, UPI000B49FC21, UPI000463D1A0, UPI000CEC0130, UPI000A018989, 
#  A0A2S8DNM3, A0A366VEF5, A0A0N5C5A0, A0A0N4Z2Y8, UPI00097154A5, A0A0P5BT45, A0A2X1KLB6, UPI00024948C4, T1WBY8, 
#  A0A0P5IY05, UPI000C16032A, UPI000C15FEF7, A0A1T8W564, A0A3F3DV56, UPI0009A63133, UPI000D1E757E, A0A0N5C599, 
#  UPI000B49D86F, UPI00093114A1, UPI000987582D, UPI000D52945D, S7JXP7, A0A0P6AAT9, UPI000E1CC48A, UPI00076B9105, 
#  A0A3A0MFJ5, L2T4X8, D1KAJ3, A0A0P5VXD9, UPI000B94675D, A0A181ZKH8, UPI0009B7AE35, UPI0005CE8979, UPI00067DD794, 
#  A0A2H1T2D3, UPI000515E2B2, UPI000935EE41

mapping$To <- mapping$To %>%
  str_split("_") %>%
  sapply(`[[`, 1) %>%
  toupper()
write.csv(mapping,paste(exportDir,"/mapping.csv",sep=""),row.names=F,quote=F)

mappedabundanceMat <- abundanceMat %>%
  rownames_to_column("Uniprot") %>%
  left_join(mapping, by = c("Uniprot" = "From")) %>%
  dplyr::select(-Uniprot) %>%
  group_by(To) %>%
  summarize_all(sum) %>%
  ungroup() %>%
  drop_na() %>%
  column_to_rownames("To")


# filter go terms
GOTermsFiltered <- GOTerms
# Keeps only genes that has at least 1 gene in the dataset
GOTermsFiltered$genesets <- GOTermsFiltered$genesets[sapply(GOTermsFiltered$genesets, function(geneset) any(geneset %in% rownames(mappedabundanceMat)))]
# keep only GO terms with at least 1 genes
GOTermsFiltered$genesets <- GOTermsFiltered$genesets[sapply(GOTermsFiltered$genesets, length) >= 5]
GOTermsFiltered$names <- GOTermsFiltered$names[names(GOTermsFiltered$genesets)]

# This is specific to the GLDS-212 dataset
# harmonize the sample name
colnames(mappedabundanceMat) <- colnames(mappedabundanceMat) %>%
  str_replace("GLDS-212_metagenomics_", "") %>%
  str_replace("_S[0-9]+_L001_R1_001_seqNames_kneaddata_Abundance.RPKs", "") %>%
  make.names() %>%
  toupper()
write.csv(mappedabundanceMat,paste(exportDir,"/mappedAbundanceMatrix.csv",sep=""),row.names=T,quote=F)

conditions = list(
  "basal-control" = c("NASA1.M1.FEC1", "NASA2.M2.FEC2", "NASA3.M3.FEC2", "NASA4.M4.FEC2", "NASA5.M5.FEC1", "NASA6.M6.FEC3", "NASA7.M7.FEC2", "NASA8.M8.FEC2", "NASA9.M9.FEC1", "NASA10.M10.FEC2", "NASA11.M12.FEC3"),
  "vivarium-control" = c("NASA12.M13.FEC1", "NASA13.M14.FEC3", "NASA14.M15.FEC1", "NASA15.M16.FEC2", "NASA16.M18.FEC3", "NASA17.M19.FEC2", "NASA18.M20.FEC3"),
  "spaceflight" = c("NASA19.M23.FEC2", "NASA20.M24.FEC3", "NASA21.M25.FEC2", "NASA22.M26.FEC1", "NASA32.M28.FEC4", "NASA23.M29.FEC3", "NASA24.M30.FEC3")
)

# run DE analysis
p_load(SummarizedExperiment)

comparisons = expand.grid(
  control = c("basal-control", "vivarium-control"),
  spaceflight = "spaceflight",
  stringsAsFactors = FALSE
)

allDEResults <- lapply(seq_len(nrow(comparisons)), function(i) {
  controlSamples <- conditions[[comparisons$control[i]]]
  spaceflightSamples <- conditions[[comparisons$spaceflight[i]]]

  expr <- cbind(mappedabundanceMat[, controlSamples], mappedabundanceMat[, spaceflightSamples])
  # expr <- t(t(expr) / colSums(expr) * 1e6)
  expr <- log2(expr + 1) %>% as.matrix()

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = data.frame(
      condition = rep(c("control", "spaceflight"), c(length(controlSamples), length(spaceflightSamples))),
      row.names = colnames(expr)
    )
  )

  design <- model.matrix(~0 + condition, data = SummarizedExperiment::colData(se))
  contrast <- limma::makeContrasts(conditionspaceflight - conditioncontrol, levels = design)

  result <- RCPA::runDEAnalysis(
    se,
    method = "limma",
    design = design,
    contrast = contrast,
    annotation = data.frame(FROM = rownames(se), TO = rownames(se), stringsAsFactors = FALSE)
  )
## Warning: Zero sample variances detected, have been offset away from zero
#  Likely due to zeros in both treatments - do we need to filter before running? What impact would that have downstream?  
  
  list(
    comparison = paste(comparisons$spaceflight[i], comparisons$control[i], sep = "_vs_"),
    result = result
  )
})

for (deRes in allDEResults){
  res <- rowData(deRes$result) %>% as.data.frame()
  write.csv(file = file.path("exportKnead", paste0(deRes$comparison, ".csv")), res, row.names = TRUE)
}

# run geneset analysis
methodsToRun <- c("ora", "fgsea", "ks", "wilcox")

allGSResults <- lapply(allDEResults, function(DEResult) {
  lapply(methodsToRun, function(method) {
    result <- RCPA::runGeneSetAnalysis(
      DEResult$result,
      genesets = GOTermsFiltered,
      method = method,
      FgseaArgs = list(sampleSize = nrow(SummarizedExperiment::colData(DEResult$result)))
    )

    result$comparison <- DEResult$comparison
    result$method <- method

    result
  }) %>% bind_rows()

}) %>% bind_rows()

# Why commented? write.csv(file = file.path("exportKnead","goterms-results.csv"), as.data.frame(allGSResults), row.names = TRUE)
write.csv(file = file.path(exportDir,"goterms-results.csv"), as.data.frame(allGSResults), row.names = TRUE)

consensusResult <- allGSResults %>% group_by(ID) %>% summarize(
  p.fisher = RCPA:::.runFisher(p.value),
  p.stouffer = RCPA:::.runStouffer(p.value),
) %>% left_join(allGSResults, by = "ID") %>% dplyr::select(ID, name, p.fisher) %>% as.data.frame() %>% unique()

# Why commented? write.csv(file = file.path("exportKnead","consensus-results.csv"), as.data.frame(consensusResult), row.names = TRUE)
write.csv(file = file.path(exportDir,"consensus-results.csv"), as.data.frame(consensusResult), row.names = TRUE)

## Added from 5_Microbiome-Analysis.R
metaRes <- RCPA::runPathwayMetaAnalysis(
  allGSResults %>%
    group_by(comparison) %>%
    group_split(),
  method = "REML"
)

metaRes$name <- metaRes$name$name
metaRes$pathwaySize <- metaRes$pathwaySize$pathwaySize

write.csv(file = file.path(exportDir,"meta-results.csv"), as.data.frame(metaRes), row.names = TRUE)


