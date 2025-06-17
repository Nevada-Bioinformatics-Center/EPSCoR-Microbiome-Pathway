#!/usr/bin/env Rscript
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#   Generate GO Terms - for downstream CPA analysis          #
#   Authors: Kanishka Manna                                  #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
goterms_file <- args[1]
gene2go_file <- args[2]
geneinfo_file <- args[3]
goobo_file <- args[4]


# Install and load required packages
library(pacman)
p_load(tidyverse, RCurl, data.table, stringr)


# Read gene2go
cat("Reading gene2go.gz ...\n")
gene2go <- fread(gene2go_file, sep = "\t", header = TRUE, fill = TRUE, showProgress = FALSE)

# Use 'Category' column for namespace
cat("Filtering for biological process ...\n")
gene2go <- gene2go[Category == "Process"]

# Read gene_info
cat("Reading All_Data.gene_info.gz ...\n")
geneinfo <- fread(geneinfo_file, select = c("GeneID", "Symbol"), showProgress = FALSE)


# Parse go.obo for GO term names
cat("Parsing go.obo ...\n")
obo_lines <- readLines(goobo_file)
term_blocks <- strsplit(paste(obo_lines, collapse = "\n"), "\n\\[Term\\]\n")[[1]]
term_blocks <- term_blocks[grepl("^id: GO:", term_blocks)]
go_names <- sapply(term_blocks, function(block) {
  id <- str_match(block, "id: (GO:\\d+)")[,2]
  name <- str_match(block, "name: ([^\n]+)")[,2]
  c(id = id, name = name)
})
go_names <- setNames(go_names["name", ], go_names["id", ])


# Map GO IDs to gene symbols
cat("Mapping GO IDs to gene symbols ...\n")
genesets <- split(gene2go$GeneID, gene2go$GO_ID)
genesets <- lapply(genesets, function(ids) {
  syms <- geneinfo[GeneID %in% ids, unique(Symbol)]
  syms <- toupper(syms)
  syms[syms != ""]
})
# Filter out small gene sets
genesets <- genesets[sapply(genesets, length) >= 5]


# Save GOTerms.rds
cat("Saving GOTerms.rds ...\n")
GOTerms <- list(
  database = "GO",
  genesets = genesets,
  names = go_names[names(genesets)]
)
saveRDS(GOTerms, goterms_file)
cat(sprintf("Done! Saved to %s\n", goterms_file))