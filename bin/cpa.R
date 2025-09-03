#!/usr/bin/env Rscript
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#   Consensus Pathway Analysis                               #
#   Authors: Melanie Hess, Cassandra K. Hui, Kanishka Manna  #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# NOTE from Kanishka
# RCPA package is not insalling properly with GitHub for Cassandra. 
# She used install.packages() to install it.
# It is installing on my workstation with devtools::install_github()
# Should take a better look into this.


# ---------- Parse command-line arguments ---------- #
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: cpa_analysis.R <samplesheet> <humann_dir> <goterms.rds>")
}
samplesheet <- args[1]
humann_dir <- args[2]
goterms_file <- args[3]

output_dir <- "cpa_out"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- Load libraries ---------- #
options(repos = c(CRAN = "https://cloud.r-project.org"))
pacman::p_load(tidyverse, RCurl, R.utils, fgsea, pbmcapply, stringr, devtools, BiocManager, 
               UniProt.ws, stringdist, SummarizedExperiment, gridpattern, ggpattern, units, sf, igraph, 
               ggraph, RColorBrewer, ggplot2)
               
devtools::install_github("tinnlab/RCPA")
library(RCPA)


# ---------- Load metadata ---------- #
meta <- read.csv(samplesheet, sep = ",", stringsAsFactors = FALSE)
if (!"exp_conditions" %in% colnames(meta)) stop("No 'exp_conditions' column found in samplesheet.")
factor <- meta$exp_conditions


# ---------- Load GO Terms ---------- #
GOTerms <- readRDS(goterms_file)


# ---------- Load HUMAnN genefamilies abundance files ---------- #
gaFiles <- list.files(humann_dir, pattern = "*_genefamilies.tsv", full.names = TRUE, recursive = TRUE)
abundances <- lapply(gaFiles, function(file) {
  df <- read_tsv(file, comment = "#", col_names = FALSE, show_col_types = FALSE)
  colnames(df) <- c("Uniprot", gsub("genefamilies.tsv", "Abundance-RPKs", basename(file)))
  df
})
allGenes <- Reduce(union, lapply(abundances, `[[`, "Uniprot"))
abundanceMat <- lapply(abundances, function(df) {
  data.frame(Uniprot = allGenes) %>%
    left_join(df, by = "Uniprot") %>%
    mutate(across(-Uniprot, ~replace_na(., 0))) %>%
    dplyr::select(-Uniprot)
}) %>% bind_cols()
rownames(abundanceMat) <- allGenes
abundanceMat <- abundanceMat[!str_detect(rownames(abundanceMat), "unclassified|UNMAPPED"),]
rownames(abundanceMat) <- rownames(abundanceMat) %>%
  str_replace("^UniRef(50|90)(_ec_filtered)?_", "") %>%
  str_replace("^_", "")

write.csv(abundanceMat, file.path(output_dir, "abundances.csv"), row.names = TRUE, quote = FALSE)


# ---------- UniProt ID to Gene Name Mapping ---------- #
mapping_file <- file.path(output_dir, "mapping.csv")
if (file.exists(mapping_file)) {
  mapping <- read.csv(mapping_file, stringsAsFactors = FALSE)
  existing_mappings <- unique(mapping$From)
} else {
  mapping <- data.frame(From=character(), To=character(), stringsAsFactors = FALSE)
  existing_mappings <- character(0)
}
all_ids <- rownames(abundanceMat)
remaining_ids <- setdiff(all_ids, existing_mappings)
if (length(remaining_ids) > 0) {
  message(sprintf("%d IDs already mapped, %d remaining to process", length(existing_mappings), length(remaining_ids)))
  batch_size <- 100
  n_batches <- ceiling(length(remaining_ids) / batch_size)
  for (i in 1:n_batches) {
    start_idx <- (i-1) * batch_size + 1
    end_idx <- min(i * batch_size, length(remaining_ids))
    batch_ids <- remaining_ids[start_idx:end_idx]
    message(sprintf("Processing batch %d of %d (%d genes)", i, n_batches, length(batch_ids)))
    batch_result <- tryCatch({
      R.utils::withTimeout({
        UniProt.ws::mapUniProt(
          from = "UniProtKB_AC-ID",
          to = "Gene_Name",
          columns = character(0L),
          query = batch_ids,
          verbose = FALSE,
          debug = FALSE,
          paginate = FALSE
        )
      }, timeout = 600)
    }, error = function(e) {
      message("Error in batch ", i, ": ", e$message)
      data.frame(From = batch_ids, To = batch_ids)
    })
    mapping <- rbind(mapping, batch_result)
    mapping <- mapping[!duplicated(mapping$From), ]
    write.csv(mapping, mapping_file, row.names = FALSE, quote = FALSE)
    Sys.sleep(0.5)
  }
} else {
  message("All IDs already mapped!")
}
mapping$To <- mapping$To %>%
  str_split("_") %>%
  sapply(`[[`, 1) %>%
  toupper()
write.csv(mapping, file.path(output_dir, "mapping.csv"), row.names = FALSE, quote = FALSE)

mappedabundanceMat <- abundanceMat %>%
  rownames_to_column("Uniprot") %>%
  left_join(mapping, by = c("Uniprot" = "From")) %>%
  dplyr::select(-Uniprot) %>%
  group_by(To) %>%
  summarize_all(sum) %>%
  ungroup() %>%
  drop_na() %>%
  column_to_rownames("To")


# ---------- Filter GO Terms ---------- #
GOTermsFiltered <- GOTerms
GOTermsFiltered$genesets <- GOTermsFiltered$genesets[sapply(GOTermsFiltered$genesets, function(geneset) any(geneset %in% rownames(mappedabundanceMat)))]
GOTermsFiltered$genesets <- GOTermsFiltered$genesets[sapply(GOTermsFiltered$genesets, length) >= 5]
GOTermsFiltered$names <- GOTermsFiltered$names[names(GOTermsFiltered$genesets)]


# ---------- Harmonize sample names ---------- #
sample_names <- meta$sample
clean_column_names <- function(colnames) {
  colnames <- gsub("_Abundance.*$", "", colnames)
  colnames <- gsub("-Abundance.*$", "", colnames)
  return(colnames)
}
colnames(mappedabundanceMat) <- clean_column_names(colnames(mappedabundanceMat))

write.csv(mappedabundanceMat, file.path(output_dir, "mappedAbundanceMatrix.csv"), row.names = TRUE, quote = FALSE)


# ---------- Differential Expression Analysis ---------- #
groups <- unique(factor)
comparisons <- t(combn(groups, 2))
colnames(comparisons) <- c("group1", "group2")
comparisons <- as.data.frame(comparisons, stringsAsFactors = FALSE)

allDEResults <- lapply(seq_len(nrow(comparisons)), function(i) {
  group1 <- comparisons$group1[i]
  group2 <- comparisons$group2[i]
  group1_samples <- meta$sample[factor == group1]
  group2_samples <- meta$sample[factor == group2]
  group1_samples <- intersect(group1_samples, colnames(mappedabundanceMat))
  group2_samples <- intersect(group2_samples, colnames(mappedabundanceMat))
  expr <- cbind(mappedabundanceMat[, group1_samples, drop=FALSE],
                mappedabundanceMat[, group2_samples, drop=FALSE])
  expr <- log2(expr + 1) %>% as.matrix()

  # Variance check: skip only this comparison if all genes have zero variance
  if (all(apply(expr, 1, var) == 0)) {
    warning(sprintf("Skipping comparison %s vs %s: no genes with non-zero variance.", group2, group1))
    return(NULL)
  }

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = data.frame(
      condition = rep(c("group1", "group2"), c(length(group1_samples), length(group2_samples))),
      row.names = colnames(expr)
    )
  )
  design <- model.matrix(~0 + condition, data = SummarizedExperiment::colData(se))
  contrast <- limma::makeContrasts(conditiongroup2 - conditiongroup1, levels = design)
  result <- RCPA::runDEAnalysis(
    se,
    method = "limma",
    design = design,
    contrast = contrast,
    annotation = data.frame(FROM = rownames(se), TO = rownames(se), stringsAsFactors = FALSE)
  )
  list(
    comparison = paste(group2, "vs", group1, sep = "_"),
    result = result
  )
})

allDEResults <- allDEResults[!sapply(allDEResults, is.null)]

dir.create(file.path(output_dir, "exportKnead"), showWarnings = FALSE)
for (deRes in allDEResults){
  res <- rowData(deRes$result) %>% as.data.frame()
  write.csv(file = file.path(output_dir, "exportKnead", paste0(deRes$comparison, ".csv")), res, row.names = TRUE)
}


# ---------- Gene Set Enrichment Analysis ---------- #
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
write.csv(file = file.path(output_dir, "goterms-results.csv"), as.data.frame(allGSResults), row.names = TRUE)


# ---------- Consensus Pathway Analysis ---------- #
consensusResult <- allGSResults %>% group_by(ID) %>% summarize(
  p.fisher = RCPA:::.runFisher(p.value),
  p.stouffer = RCPA:::.runStouffer(p.value),
) %>% left_join(allGSResults, by = "ID") %>% dplyr::select(ID, name, p.fisher) %>% as.data.frame() %>% unique()
write.csv(file = file.path(output_dir, "consensus-results.csv"), as.data.frame(consensusResult), row.names = TRUE)


###############################################
# ---------- Pathway Meta-Analysis ---------- #
#n_comparisons <- length(unique(allGSResults$comparison))
#if (n_comparisons >= 2) {
#  metaRes <- RCPA::runPathwayMetaAnalysis(
#    allGSResults %>%
#      group_by(comparison) %>%
#      group_split(),
#    method = "REML"
#  )
#  metaRes$name <- metaRes$name$name
#  metaRes$pathwaySize <- metaRes$pathwaySize$pathwaySize
#  write.csv(file = file.path(output_dir, "meta-results.csv"), as.data.frame(metaRes), row.names = TRUE)
#} else {
#  warning("Meta-analysis requires results from two or more comparisons. Skipping meta-analysis step.")
#}
#############################################


# ---------- Static Plots ---------- #

# Add comparison column. Take the comparison with the lowest p.value for each ID
consensusResult <- consensusResult %>%
  left_join(
    allGSResults %>%
      group_by(ID) %>%
      slice_min(order_by = p.value, n = 1) %>%
      dplyr::select(ID, comparison),
    by = "ID"
  )

# Select top 25 GO terms by lowest p.fisher value
top_terms <- consensusResult %>%
  arrange(p.fisher) %>%
  head(25)

# Build the GO term network graph
go_graph <- graph.empty(directed = FALSE)
go_graph <- add_vertices(go_graph, nrow(top_terms),
                         name = top_terms$name,
                         category = top_terms$category,
                         significance = -log10(top_terms$p.fisher),
                         comparison = top_terms$comparison)

# Add edges between nodes sharing common words (excluding stop words)
for(i in 1:(nrow(top_terms)-1)) {
  for(j in (i+1):nrow(top_terms)) {
    # Only connect nodes from the same comparison
    if(top_terms$comparison[i] == top_terms$comparison[j]) {
      words_i <- unlist(strsplit(tolower(top_terms$name[i]), " "))
      words_j <- unlist(strsplit(tolower(top_terms$name[j]), " "))
      common_words <- intersect(words_i, words_j)
      common_words <- common_words[!common_words %in% 
                                     c("of", "the", "a", "in", "to", "by", "and", "or")]
      if(length(common_words) > 0) {
        go_graph <- add_edges(go_graph, c(i, j), weight = length(common_words), 
                              comparison = top_terms$comparison[i])
      }
    }
  }
}

# Cluster nodes using Louvain method for community detection
comm <- cluster_louvain(go_graph)
membership <- membership(comm)
# Set node and edge attributes for plotting
V(go_graph)$size <- V(go_graph)$significance * 3
V(go_graph)$community <- membership
E(go_graph)$width <- E(go_graph)$weight
# For reproducible layout
set.seed(42)
  
# Create and save the network plot
# Need a color palette with more colors or use default
network_plot <- ggraph(go_graph, layout = "fr") + 
    geom_edge_link(alpha = 0.3) +
    geom_node_point(aes(size = significance, color = as.factor(community)), alpha = 0.8) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3.5) +
    #scale_color_brewer(palette = "Set1") +
    scale_size_continuous(range = c(3, 10)) +
    facet_wrap(~ comparison) +
    labs(title = "Network of Top 25 GO Terms",
         #subtitle = "Grouped by Comparison",
         size = "-log10(p-value)",
         color = "Category") +
    theme_void() +
    theme(
      text = element_text(size = 18),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 16.5, hjust = 0.5),
      legend.position = "right",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )

ggsave(network_plot, file = file.path(output_dir,"Pathway_Networks_Plot.png"), width = 10.5, height = 7, dpi = 600)


# Create and save the bubble plot for the top GO terms
bubble_plot <- ggplot(top_terms, aes(x = comparison, y = reorder(name, -log10(p.fisher)))) +
  # geom_segment(aes(x = 0, xend = -log10(p.fisher), 
  #                  y = reorder(name, -log10(p.fisher)), 
  #                  yend = reorder(name, -log10(p.fisher))), color = "grey50") +
  geom_point(aes(size = -log10(p.fisher), color = comparison), alpha = 0.8) +
  #scale_color_viridis_c() +
  scale_color_brewer(palette = "Set1") +
  scale_size_continuous(range = c(3, 8)) +
  theme_minimal() +
  labs(title = "Top 25 Significant GO Terms",
       x = "Comparison", y = NULL,
       size = "-log10(p-value)", 
       color = "Comparison"
  ) +
  theme(
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 0.8),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 14 * 1.2),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(bubble_plot, file = file.path(output_dir,"Pathway_Bubble_Plot.png"), width = 11, height = 10, dpi = 600)



###################################################

message("Consensus Pathway Analysis complete. Results saved to: ", output_dir)

############
# THE END! #
############