#!/usr/bin/env Rscript
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#   Descriptive profiling                                    #
#   Authors: Kanishka Manna                                  #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript desc_prof_abun.R <joined_pathabundance.tsv> <top_genus_out.tsv> <top_species_out.tsv>")
}

input_file <- args[1]
genus_out <- args[2]
species_out <- args[3]

suppressPackageStartupMessages({
  library(mia)
  library(dplyr)
})

# Read the file as a plain table to check for taxonomy
tab <- read.table(input_file, header=TRUE, sep="\t", comment.char = "", check.names = FALSE)
# Remove leading '#' from column name if present
colnames(tab)[1] <- sub("^#\\s*", "", colnames(tab)[1])

# If only UNMAPPED/UNINTEGRATED, or no taxonomy in first column, skip
if (all(tab[[1]] %in% c("UNMAPPED", "UNINTEGRATED")) ||
    !any(grepl("\\|", tab[[1]]))) {
    message("No taxonomy data found. Skipping descriptive profiling.")
    file.create(genus_out)
    file.create(species_out)
    quit(save = "no", status = 0)
}

# Import HUMAnN joined pathway abundance table
se <- importHUMAnN(input_file)
abund <- assay(se)
df <- as.data.frame(rowData(se))
df <- cbind(df, as.matrix(abund))

# Ensure relevant columns are character
if ("Pathway" %in% colnames(df)) df$Pathway <- as.character(df$Pathway)
if ("genus" %in% colnames(df)) df$genus <- as.character(df$genus)
if ("species" %in% colnames(df)) df$species <- as.character(df$species)

# For each pathway, find the genus with the highest total abundance
top_genus_per_pathway <- df %>% 
  group_by(Pathway, genus) %>% 
  summarise(total_abundance = sum(across(where(is.numeric)), na.rm = TRUE), .groups = "drop") %>% 
  arrange(Pathway, desc(total_abundance)) %>% 
  group_by(Pathway) %>% 
  slice_head(n = 1) %>% ungroup() %>% 
  filter(!is.na(Pathway) & !is.na(genus) & !is.na(total_abundance))

# For each pathway, find the species with the highest total abundance
top_species_per_pathway <- df %>%
  group_by(Pathway, species) %>%
  summarise(total_abundance = sum(across(where(is.numeric)), na.rm = TRUE), .groups = "drop") %>%
  arrange(Pathway, desc(total_abundance)) %>%
  group_by(Pathway) %>%
  slice_head(n = 1) %>%
  ungroup() %>% 
  filter(!is.na(Pathway) & !is.na(species) & !is.na(total_abundance))

# Saving the descriptive result
write.table(top_genus_per_pathway, genus_out, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(top_species_per_pathway, species_out, sep = "\t", row.names = FALSE, quote = FALSE)