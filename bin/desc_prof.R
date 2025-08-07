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

# ---------- Load libraries ---------- #
options(repos = c(CRAN = "https://cloud.r-project.org"))
library(pacman)
pacman::p_load(tidyverse)


# ---------- Load merged HUMAnN pathway abundance file ---------- #
df <- read.table(input_file, header=TRUE, sep="\t", comment.char = "", check.names = FALSE)
colnames(df)[1] <- sub("^#\\s*", "", colnames(df)[1])


# If only UNMAPPED/UNINTEGRATED, or no taxonomy in first column, skip
if (all(df[[1]] %in% c("UNMAPPED", "UNINTEGRATED")) ||
    !any(grepl("\\|", df[[1]]))) {
    msg <- "No taxonomy data found. Skipping descriptive profiling."
    writeLines(msg, output_file)
    quit(save = "no", status = 0)
}


# Split the first column into taxonomy
df$taxonomy <- ifelse(grepl("\\|", df[[1]]), sub("^[^|]*\\|", "", df[[1]]), NA)

df$genus <- ifelse(!is.na(df$taxonomy) & grepl("g__", df$taxonomy),
                   sub("(g__[^.]+).*", "\\1", df$taxonomy), NA)

df$species <- ifelse(!is.na(df$taxonomy) & grepl("s__", df$taxonomy),
                     sub(".*(s__[^ ]+).*", "\\1", df$taxonomy), NA)

df$Pathway <- ifelse(grepl("\\|", df[[1]]), sub("\\|.*", "", df[[1]]), df[[1]])

# Ensure relevant columns are character
df$Pathway <- as.character(df$Pathway)
df$genus <- as.character(df$genus)
df$species <- as.character(df$species)

# Remove rows with missing or malformed taxonomy
df_clean <- df %>%
  filter(!is.na(Pathway)) %>%
  filter(!Pathway %in% c("UNMAPPED", "UNINTEGRATED")) %>%
  filter(!grepl("^\\d+(\\.\\d+)?$", Pathway)) %>%
  filter(!is.na(genus))


# Top genus per pathway
top_genus <- df_clean %>%
  filter(!is.na(genus)) %>%
  group_by(Pathway, genus) %>%
  filter(rowSums(across(where(is.numeric)), na.rm = TRUE) > 0) %>%
  summarise(GenusCPM = max(across(where(is.numeric)), na.rm = TRUE), .groups = "drop") %>%
  arrange(Pathway, desc(GenusCPM)) %>%
  group_by(Pathway) %>%
  slice_head(n = 1) %>%
  ungroup()

# Top species per pathway
#top_species <- df %>%
#  filter(!is.na(species)) %>%
#  group_by(Pathway, species) %>%
#  summarise(SpeciesCPM = sum(across(where(is.numeric)), na.rm = TRUE), .groups = "drop") %>%
#  arrange(Pathway, desc(SpeciesCPM)) %>%
#  group_by(Pathway) %>%
#  slice_head(n = 1) %>%
#  ungroup()

# Merge top genus and species for each pathway
#final <- full_join(top_genus, top_species, by = "Pathway") %>%
#  mutate(
#    CountsPerMillion = coalesce(SpeciesCPM, GenusCPM)
#  ) %>%
#  select(Pathway, genus, species, CountsPerMillion)

colnames(top_genus) <- c("Pathway", "Genus", "CountsPerMillion")
write.table(top_genus, output_file, sep = "\t", row.names = FALSE, quote = FALSE)

# THE END!
##########