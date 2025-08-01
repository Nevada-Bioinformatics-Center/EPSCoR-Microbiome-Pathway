#!/usr/bin/env Rscript
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#   Descriptive profiling                                    #
#   Authors: Kanishka Manna                                  #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript desc_prof_abun.R <joined_pathabundance.tsv> <output.tsv>")
}

input_file <- args[1]
output_file <- args[2]

suppressPackageStartupMessages({
  library(mia)
  library(dplyr)
})

# Read the file as a plain table to check for taxonomy
tab <- read.table(input_file, header=TRUE, sep="\t", comment.char = "", check.names = FALSE)
colnames(tab)[1] <- sub("^#\\s*", "", colnames(tab)[1])

# If only UNMAPPED/UNINTEGRATED, or no taxonomy in first column, skip
if (all(tab[[1]] %in% c("UNMAPPED", "UNINTEGRATED")) ||
    !any(grepl("\\|", tab[[1]]))) {
    msg <- "No taxonomy data found. Skipping descriptive profiling."
    writeLines(msg, output_file)
    quit(save = "no", status = 0)
}

# Import HUMAnN joined pathway abundance table
se <- importHUMAnN(input_file)
abund <- assay(se)
df <- as.data.frame(rowData(se))
df <- cbind(df, as.matrix(abund))

# Split the first column into pathway and taxonomy
df$Pathway <- sub("\\|.*", "", df[[1]])
df$taxonomy <- sub("^[^|]*\\|?", "", df[[1]])

# Extract genus and species from taxonomy
df$genus <- sub("(g__[^.]+).*", "\\1", df$taxonomy)
df$species <- sub(".*(s__[^ ]+).*", "\\1", df$taxonomy)
df$genus[!grepl("g__", df$genus)] <- NA
df$species[!grepl("s__", df$species)] <- NA

# Ensure relevant columns are character
df$Pathway <- as.character(df$Pathway)
df$genus <- as.character(df$genus)
df$species <- as.character(df$species)

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

# Combine genus and species results for each pathway
combined <- merge(
  top_genus_per_pathway[, c("Pathway", "genus", "total_abundance")],
  top_species_per_pathway[, c("Pathway", "species", "total_abundance")],
  by = "Pathway",
  suffixes = c("_genus", "_species")
)

# Choose the higher abundance between genus and species for each pathway
combined$abundance <- pmax(combined$total_abundance, combined$total_abundance_species, na.rm = TRUE)

# Select only relevant columns
final <- combined[, c("Pathway", "genus", "species", "abundance")]

# Save to a single output file
write.table(final, output_file, sep = "\t", row.names = FALSE, quote = FALSE)