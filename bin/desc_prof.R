#!/usr/bin/env Rscript
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#   Descriptive profiling                                    #
#   Authors: Kanishka Manna, Hans Vasquez-Gross &            #
#            Cassandra K. Hui                                #
#   Updated: fix list-column bug; use per-sample MAX metric  #
#   Updated: Added exp conditions, sample columns; fix       #
#            abundance values from unclassified and          #
#            classified taxa                                 #

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript desc_prof.R <joined_pathabundance.tsv> <output1.tsv> <output2.tsv> <output3.tsv> <output4.tsv> <samplesheet.csv>")
}

input_file  <- args[1]
output_file1 <- args[2]
output_file2 <- args[3]
output_file3 <- args[4]
output_file4 <- args[5]
samplesheet <- args[6]


# input_file <- "joined_norm_pathabundance.tsv"
# output_file1 <- "pathways.tsv"
# output_file2 <- "sample_statistics.tsv"
# output_file3 <- "pathways_by_genus.tsv"
# output_file4 <- "genus_statistics.tsv"
# 
# samplesheet <- "seedfile.csv"

# ---------- Libraries ----------
options(repos = c(CRAN = "https://cloud.r-project.org"))
suppressPackageStartupMessages({
  library(pacman)
})
pacman::p_load(data.table, dplyr, tidyr, stringr)

# ---------- Load sample metadata ----------
meta <- data.table::fread(samplesheet, sep = ",", header = TRUE, data.table = FALSE)
meta <- meta[, intersect(c("sample", "exp_conditions"), colnames(meta)), drop = FALSE]

# ---------- Load the merged renorm path-abundance file ----------
df <- data.table::fread(
  input_file,
  sep = "\t",
  header = TRUE,
  quote = "",
  fill = TRUE,
  data.table = FALSE,
  showProgress = FALSE
)

# Normalize the first header (some HUMAnN files start with "#")
colnames(df)[1] <- sub("^#\\s*", "", colnames(df)[1])


##############################
# Separate the two dataframes
##############################

# Extract total (unstratified) pathways
df_total_pathways <- df %>% filter(!grepl("\\|", .[[1]]))

# Extract stratified (taxa-associated) pathways
df_stratified <- df %>% filter(grepl("\\|", .[[1]]))

# ---------- Parse taxonomy + pathway ----------

df_stratified <- df_stratified %>%
  mutate(
    taxonomy = ifelse(grepl("\\|", .[[1]]), sub("^[^|]*\\|", "", .[[1]]), NA),
    # genus    = grepl("g__", taxonomy),
    #                   sub("(g__[^| ]+).*", "\\1", taxonomy),
    genus    = sub("\\.s__.*", "", taxonomy),  # <-- This removes the species/strain part
    Pathway  = ifelse(grepl("\\|", .[[1]]), sub("\\|.*", "", .[[1]]), .[[1]])
  ) %>%
  mutate(
    Pathway = as.character(Pathway),
    taxonomy = as.character(taxonomy),
    genus   = as.character(genus)
  )


## Sum rows and sort by size

df_sum <- df_total_pathways %>%
  mutate(Total = rowSums(select_if(., is.numeric), na.rm = TRUE)) %>%
  arrange(desc(Total))


df_sum_stratified <- df_stratified %>%
  mutate(Total = rowSums(select_if(., is.numeric), na.rm = TRUE)) %>%
  arrange(desc(Total))

# ---------- Create Sample Information Dataframe ----------

# Identify sample columns (before any filtering)
annot_cols_initial <- c("taxonomy", "genus", "Pathway")
sample_cols_initial <- setdiff(names(df_sum), annot_cols_initial)

# Calculate total abundance per sample
total_abundance_per_sample <- colSums(df_sum[sample_cols_initial], na.rm = TRUE)

# Calculate UNMAPPED/UNINTEGRATED abundance per sample
unmapped_df <- df_sum %>% filter(Pathway %in% "UNMAPPED")
unmapped_abundance_per_sample <- colSums(unmapped_df[sample_cols_initial], na.rm = TRUE)


unintegrated_df <- df_sum %>% filter(Pathway %in% "UNINTEGRATED")
unintegrated_abundance_per_sample <- colSums(unintegrated_df[sample_cols_initial], na.rm = TRUE)

# Create sample information dataframe
sample_info <- data.frame(
  Sample = sample_cols_initial,
  Total_Abundance_CPM = total_abundance_per_sample,
  Unmapped_CPM = unmapped_abundance_per_sample,
  Unmapped_Percent =
    round((unmapped_abundance_per_sample / total_abundance_per_sample) * 100, 2),
  Unintegrated_CPM = unintegrated_abundance_per_sample,
  Unintegrated__Percent =
    round((unintegrated_abundance_per_sample / total_abundance_per_sample) * 100, 2)
)

# Drop non-informative rows after keeping data
df_clean <- df_sum %>%
  filter(
    !is.na(Pathway),
    !Pathway %in% c("UNMAPPED", "UNINTEGRATED")
  )

df_clean_stratified <- df_sum_stratified %>%
  filter(
    !is.na(Pathway),
    !Pathway %in% c("UNMAPPED", "UNINTEGRATED"),
    !taxonomy %in% "unclassified"
  )


# ---------- Add Pathway Statistics Per Sample to Sample Info ----------

# Calculate pathway statistics per sample
pathway_stats_per_sample <- data.frame(
  Sample = sample_cols_initial,
  Pathways_Detected = sapply(sample_cols_initial, function(col) {
    sum(df_clean[[col]] > 0, na.rm = TRUE)
  }),
  Pathways_with_Genus = sapply(sample_cols_initial, function(col) {
    # Get pathways with abundance > 0 in this sample
    active_pathways <- df_clean_stratified[[col]] > 0
    # Count unique pathways that are active
    length(unique(df_clean_stratified$Pathway[active_pathways]))
  }),
  Unique_Genera_Detected = sapply(sample_cols_initial, function(col) {
    active_pathways <- df_clean_stratified[[col]] > 0
    length(unique(df_clean_stratified$genus[active_pathways]))
  }),
  Unique_Species_Detected = sapply(sample_cols_initial, function(col) {
    active_pathways <- df_clean_stratified[[col]] > 0
    length(unique(df_clean_stratified$taxonomy[active_pathways]))
  })
)



# Merge with existing sample_info
sample_info <- merge(sample_info, pathway_stats_per_sample, by = "Sample", all = TRUE)

# Grouping by experimental conditions
# Remove '_Abundance' suffix from Sample column before joining
sample_info <- sample_info %>%
  mutate(Sample_clean = sub("_Abundance$", "", Sample))

# Now join using the cleaned sample name
sample_info <- sample_info %>%
  left_join(meta, by = c("Sample_clean" = "sample")) %>%
  select(-Sample_clean)  # Optionally remove the helper column

# Rename the columns
colnames(sample_info) <- c(
  "Sample",
  "TotalAbundanceCPM",
  "UnmappedCPM",
  "UnmappedPercent",
  "UnintegratedCPM",
  "UnintegratedPercent",
  "PathwaysDetected",
  "PathwaysWithGenus",
  "UniqueGeneraDetected",
  "UniqueSpeciesDetected",
  "ExperimentalGroup"
)


# ---------- Prepare final table ----------
path_tbl <- df_clean %>%
  transmute(
    Pathway = Pathway,
    #Genus = taxonomy,
    Abundance = Total
  )


genus_tbl <- df_clean_stratified %>%
  transmute(
    Pathway = Pathway,
    Taxonomy = taxonomy,
    Genus = genus,
    Abundance = Total
  )

# Group by exp_conditions
# Join genus_tbl with meta to add experimental group and sample name info
genus_tbl_grouped <- genus_tbl %>%
  left_join(meta, by = character())



# # Safety guard: no list-columns allowed (keeps TSV tidy)
# stopifnot(!any(vapply(final_tbl, is.list, logical(1)))). ### What is this???


#############
# Genus here has 3565, should we cut at all? Or do we want to output full table?

#############

# Show which genera are driving which pathways
# Grouped by Genus and exp_conditions
genus_pathway_summary <- genus_tbl_grouped %>%
  group_by(Genus, exp_conditions) %>%
  summarise(
    Pathway_Count = length(unique(Pathway)),
    Total_Abundance = sum(Abundance),
    Percentage_of_Classified_Pathways = round((Total_Abundance / sum(genus_tbl$Abundance)) * 100, 2),
    Top_Pathway = Pathway[which.max(Abundance)],
    .groups = "drop"
  ) %>%
  arrange(exp_conditions, desc(Total_Abundance))

# rename columns
  colnames(genus_pathway_summary) <- c("Genus", "ExperimentalGroup", "PathwayCount",
                                     "TotalAbundance", "Percentage of Classified Pathways",
                                     "TopPathways")



#############
## Need to rename these headers too



# ---------- Write ----------
if (nrow(path_tbl) == 0L) {
  message("desc_prof: no pathways found; writing message file")
  writeLines("No pathways found", con = output_file1)
} else {
  write.table(path_tbl, output_file1, sep = "\t", row.names = FALSE, quote = FALSE)
}
if (nrow(sample_info) == 0L) {
  message("desc_prof: no sample info found; writing message file")
  writeLines("No sample info found", con = output_file2)
} else {
  write.table(sample_info, output_file2, sep = "\t", row.names = FALSE, quote = FALSE)
}
if (nrow(genus_tbl) == 0L) {
  message("desc_prof: no genus info found; writing message file")
  writeLines("No genus info found", con = output_file3)
} else {
  write.table(genus_tbl, output_file3, sep = "\t", row.names = FALSE, quote = FALSE)
}
if (nrow(genus_pathway_summary) == 0L) {
  message("desc_prof: no genus pathways found; writing message file")
  writeLines("No genus pathways found", con = output_file4)
} else {
  write.table(genus_pathway_summary, output_file4, sep = "\t", row.names = FALSE, quote = FALSE)
}
