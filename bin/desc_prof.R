#!/usr/bin/env Rscript
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#   Descriptive profiling                                    #
#   Authors: Kanishka Manna & Hans Vasquez-Gross             #
#   Updated: fix list-column bug; use per-sample MAX metric  #

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript desc_prof.R <joined_pathabundance.tsv> <output.tsv>")
}

input_file  <- args[1]
output_file <- args[2]

# ---------- Libraries ----------
options(repos = c(CRAN = "https://cloud.r-project.org"))
suppressPackageStartupMessages({
  library(pacman)
})
pacman::p_load(data.table, dplyr, tidyr, stringr)

# ---------- Load (robust) ----------
# Disable quoting; allow ragged rows to fill with NA
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

# ---------- Early exit if no taxonomy ----------
first_col <- df[[1]]
if (all(first_col %in% c("UNMAPPED", "UNINTEGRATED")) || !any(grepl("\\|", first_col))) {
  writeLines("No taxonomy data found. Skipping descriptive profiling.", output_file)
  quit(save = "no", status = 0)
}

# ---------- Parse taxonomy + pathway ----------
df <- df %>%
  mutate(
    taxonomy = ifelse(grepl("\\|", .[[1]]), sub("^[^|]*\\|", "", .[[1]]), NA),
    genus    = ifelse(!is.na(taxonomy) & grepl("g__", taxonomy),
                      sub("(g__[^| ]+).*", "\\1", taxonomy), NA),
    species  = ifelse(!is.na(taxonomy) & grepl("s__", taxonomy),
                      sub(".*(s__[^| ]+).*", "\\1", taxonomy), NA),
    Pathway  = ifelse(grepl("\\|", .[[1]]), sub("\\|.*", "", .[[1]]), .[[1]])
  ) %>%
  mutate(
    Pathway = as.character(Pathway),
    genus   = as.character(genus),
    species = as.character(species)
  )

# Drop non-informative rows
df_clean <- df %>%
  filter(
    !is.na(Pathway),
    !Pathway %in% c("UNMAPPED", "UNINTEGRATED"),
    !grepl("^\\d+(\\.\\d+)?$", Pathway),  # accidental numeric-only entries
    !is.na(genus)
  )

if (nrow(df_clean) == 0) {
  writeLines("No rows with genus-level taxonomy after filtering. Skipping descriptive profiling.", output_file)
  quit(save = "no", status = 0)
}

# ---------- Identify sample columns & coerce to numeric ----------
annot_cols  <- c("taxonomy", "genus", "species", "Pathway")
sample_cols <- setdiff(names(df_clean), annot_cols)

# Coerce candidates to numeric (non-numeric -> NA)
df_clean[sample_cols] <- lapply(df_clean[sample_cols], function(x) suppressWarnings(as.numeric(x)))

# Keep only truly numeric columns
is_num <- vapply(df_clean[sample_cols], is.numeric, logical(1))
sample_cols <- sample_cols[is_num]

if (length(sample_cols) == 0) {
  writeLines("No numeric sample columns detected. Skipping descriptive profiling.", output_file)
  quit(save = "no", status = 0)
}

# ---------- Row-wise max across samples, drop zero/NA rows ----------
row_max <- do.call(pmax, c(df_clean[sample_cols], list(na.rm = TRUE)))
df_clean$RowMax <- row_max
df_clean <- df_clean %>% filter(is.finite(RowMax), RowMax > 0)

if (nrow(df_clean) == 0) {
  writeLines("No non-zero abundances after filtering. Skipping descriptive profiling.", output_file)
  quit(save = "no", status = 0)
}

# ---------- Top genus per pathway using MAX ----------
top_genus <- df_clean %>%
  group_by(Pathway, genus) %>%
  summarise(CountsPerMillion = max(RowMax, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(CountsPerMillion)) %>%
  group_by(Pathway) %>%
  slice_max(order_by = CountsPerMillion, n = 1, with_ties = FALSE) %>%
  ungroup()

# Safety guard: no list-columns allowed (keeps TSV tidy)
stopifnot(!any(vapply(top_genus, is.list, logical(1))))

# ---------- Write ----------
colnames(top_genus) <- c("Pathway", "Genus", "CountsPerMillion")
write.table(top_genus, output_file, sep = "\t", row.names = FALSE, quote = FALSE)# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
