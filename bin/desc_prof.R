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

# Import HUMAnN joined pathway abundance table
se <- importHUMAnN(input_file)

# Get the abundance matrix and pathway-taxa info & combine 
abund <- assay(se)
df <- as.data.frame(rowData(se))
df <- cbind(df, as.matrix(abund))

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