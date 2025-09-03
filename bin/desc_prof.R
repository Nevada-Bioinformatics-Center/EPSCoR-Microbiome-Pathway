#!/usr/bin/env Rscript
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#   Descriptive profiling                                    #
#   Authors: Kanishka Manna & Hans Vasquez-Gross             #
#   Updated: fix list-column bug; use per-sample MAX metric  #
#   Updated: Added exp conditions, sample columns; fix       #
#            abundance values from unclassified and          #
#            classified taxa                                 #

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript desc_prof.R <joined_pathabundance.tsv> <output.tsv> <samplesheet.csv>")
}

input_file  <- args[1]
output_file <- args[2]
samplesheet <- args[3]

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

# ---------- Parse taxonomy + pathway ----------
df <- df %>%
  mutate(
    taxonomy = ifelse(grepl("\\|", .[[1]]), sub("^[^|]*\\|", "", .[[1]]), NA),
    genus    = ifelse(!is.na(taxonomy) & grepl("g__", taxonomy),
                      sub("(g__[^| ]+).*", "\\1", taxonomy), "unclassified"),
    Pathway  = ifelse(grepl("\\|", .[[1]]), sub("\\|.*", "", .[[1]]), .[[1]])
  ) %>%
  mutate(
    Pathway = as.character(Pathway),
    genus   = as.character(genus)
  )

# ---------- Identify sample columns & coerce to numeric ----------
annot_cols  <- c("taxonomy", "genus", "Pathway")
sample_cols <- setdiff(names(df), annot_cols)
df[sample_cols] <- lapply(df[sample_cols], function(x) suppressWarnings(as.numeric(x)))
is_num <- vapply(df[sample_cols], is.numeric, logical(1))
sample_cols <- sample_cols[is_num]

if (length(sample_cols) == 0) {
  writeLines("No numeric sample columns detected. Skipping descriptive profiling.", output_file)
  quit(save = "no", status = 0)
}

# ---------- Row-wise max across samples ----------
row_max <- do.call(pmax, c(df[sample_cols], list(na.rm = TRUE)))
max_sample_idx <- apply(df[sample_cols], 1, function(x) {
  if (all(is.na(x))) return(NA)
  which.max(x)
})
max_sample <- names(df[sample_cols])[max_sample_idx]
df$Abundance <- row_max
df$SampleName <- max_sample

# ---------- Prepare final table ----------
final_tbl <- df %>%
  transmute(
    Pathway = Pathway,
    Genus = taxonomy,
    Abundance = Abundance,
    SampleName = sub("_Abundance$", "", SampleName)
  )

# Map ExperimentalConditions from meta using SampleName, only if exp_conditions exists
if ("exp_conditions" %in% colnames(meta)) {
  final_tbl$ExperimentalConditions <- meta$exp_conditions[match(final_tbl$SampleName, meta$sample)]
}

final_tbl <- final_tbl %>%
  filter(!is.na(Genus) | Pathway == "UNMAPPED")

# Safety guard: no list-columns allowed (keeps TSV tidy)
stopifnot(!any(vapply(final_tbl, is.list, logical(1))))

# ---------- Write ----------
write.table(final_tbl, output_file, sep = "\t", row.names = FALSE, quote = FALSE)