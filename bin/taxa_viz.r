#!/usr/bin/env Rscript
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#   Visualizations of Taxa and QC                            #
#   Authors: Cassandra K. Hui, Kanishka Manna                #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
merged_taxa_file <- args[1]
samplesheet_file <- args[2]

exportDir <- "plots"
dir.create(exportDir, showWarnings = FALSE, recursive = TRUE)


# ---------- Load libraries ---------- #
options(repos = c(CRAN = "https://cloud.r-project.org"))
library(pacman)
pacman::p_load(mia, phyloseq, vegan, ggplot2, dplyr)

# ---------- Load samplesheet ---------- #
metaData <- read.csv(samplesheet_file, sep = ",", stringsAsFactors = FALSE)
rownames(metaData) <- metaData$sample

# Color by experimental condition if it exists and has more than one unique value
if ("exp_conditions" %in% colnames(metaData) && length(unique(metaData$exp_conditions)) > 1) {
  factor = "exp_conditions"
} else {
  factor = "sample"
}


# ---------- Load merged Metaphlan file ---------- #
# NOTE: Use merge_metaphlan_tables.py to generate merged_taxa_file before running this script
# Example: python merge_metaphlan_tables.py *_profile.tsv > merged_abundance_table.tsv
tse <- mia::importMetaPhlAn(merged_taxa_file, assay.type = "metaphlan")

# 'assasy.type' is set to default - "metaphlan"
# To find assay.type, run assayName(tse) in the console
phylo <- mia::convertToPhyloseq(tse, assay.type = "metaphlan")

# Add metadata to phyloseq object
sample_data(phylo) <- metaData

### Can't do a rarefaction curve without raw data (normalized doesn't work)


# ---------- Alpha Diversity Plots ---------- #

# set global theme for ggplot2
#theme_set(theme_bw())
# Estimate richness (Shannon and Inverse Simpson indices)
alpha_df <- estimate_richness(phylo, measures=c("Shannon", "InvSimpson"))
alpha_df$sample <- rownames(alpha_df)
alpha_df$exp_conditions <- as.factor(metaData$exp_conditions[match(alpha_df$sample, metaData$sample)])

# Shannon diversity plot
alpha_plot <- ggplot(alpha_df) +
  geom_point(aes(x = sample, y = Shannon, fill = .data[[factor]]), pch = 21, color = "black", size = 4, alpha = 0.8) +
  scale_fill_viridis_d(option="turbo") +
  theme_bw(base_size = 14) +
  labs(x = "Samples",
       y = "Shannon Index") +
  theme(
    axis.text.x =element_text(angle=90),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    #,
    #panel.grid = element_blank()
  )

#alpha_plot

ggsave(alpha_plot, file = file.path(exportDir, "alpha_diversity_shannon.png"), height = 6, width = 8, units = "in")


# Inverse Simpson diversity plot
alpha_plotIS <- ggplot(alpha_df) +
  geom_point(aes(x = sample, y = InvSimpson, fill = .data[[factor]]), pch = 21, color = "black", size = 4, alpha = 0.8) +
  scale_fill_viridis_d(option="turbo") +
  theme_bw(base_size = 14) +
  labs(x = "Samples",
       y = "Inverse Simpson Index") +
  theme(
    axis.text.x =element_text(angle=90),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    #,
    #panel.grid = element_blank()
  )

#alpha_plotIS

ggsave(alpha_plotIS, file = file.path(exportDir, "alpha_diversity_inverse_simpson.png"), height = 6, width = 8, units = "in")


# ---------- Top 10 Phyla Stacked Bar Plot ---------- #
phylum<-tax_glom(phylo,taxrank = "phylum")
top10<- names(sort(taxa_sums(phylum), TRUE)[1:10])
dat<-psmelt(phylum)
dat<-mutate(dat, phylum= ifelse(OTU %in% top10, paste0(phylum), "Other"))

top10plot <- ggplot(dat, aes(x=Sample,y=Abundance, col=phylum))+
  geom_col(aes( fill=phylum), color = "black", position="stack", linewidth = 0.2)+
  #theme_pubr(legend="right")+
  scale_fill_viridis_d(option="turbo") + 
  labs(x = "Samples") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x =element_text(angle=90),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.grid = element_blank()
  )
#top10plot

ggsave(top10plot, file = file.path(exportDir, "scaled_abundances_top10.png"), height = 6, width = 8, units = "in")


# ---------- NMDS Ordination Plot ---------- #
# Remove samples with all-zero or NA abundance
phylo <- prune_samples(sample_sums(phylo) > 0, phylo)
ord.nmds.bray <- ordinate(phylo, method="NMDS", distance="bray")
ordination_plot <- plot_ordination(phylo, ord.nmds.bray, color="sample", title="Bray NMDS")
ggsave(ordination_plot, file = file.path(exportDir, "Bray_NMDS_sample.png"), height = 6, width = 8, units = "in")

nmds_df <- plot_ordination(phylo, ord.nmds.bray, justDF = TRUE)

if ("exp_conditions" %in% colnames(nmds_df)) {
  nmds_df$exp_conditions <- as.factor(nmds_df$exp_conditions)
}

nmds_plot <- ggplot(nmds_df) +
  # Dashed lines for improved readability
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8, color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8, color = "gray50") +
  
  geom_point(aes(x = NMDS1, y = NMDS2, fill = .data[[factor]]), pch = 21, color = "black", size = 4, alpha = 0.8) +
  #theme_pubr(legend = "right") +
  
  # # Add ellipses with soft transparency
  # stat_ellipse(aes(x = NMDS1, y = NMDS2, col = Sample.Name, fill = Sample.Name),
  #              geom = "polygon", type = "t", alpha = 0.05,
  #              show.legend = FALSE, level = 0.95) +
  # 
  # Apply custom colors
  #scale_color_manual(values = treatment_colors) +
  scale_fill_viridis_d(option="turbo") + 
  
  # Improve overall plot aesthetics
  labs(title = "NMDS Ordination",
       caption = paste0("Stress = ", round(ord.nmds.bray$stress, 4))) +

  # Adjust theme for better readability
  theme_bw(base_size = 14) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.grid = element_blank()
  )

#nmds_plot

ggsave(nmds_plot, file = file.path(exportDir, "NMDS.png"), height = 6, width = 8, units = "in")


# ---------- Phylum Abundance Bar Plot ---------- #
phy_plot <- plot_bar(phylo, x="sample", fill="phylum") + 
  #facet_wrap(~Treatment, scales="free_x", nrow = 1) + 
  scale_fill_viridis_d(option="turbo") +
  theme_bw(base_size = 14) +
  labs(x = "Samples") +
  theme(
    axis.text.x =element_text(angle=90),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.grid = element_blank()
  )

#phy_plot

ggsave(phy_plot, file = file.path(exportDir, "phylum_abundence.png"), height = 6, width = 8, units = "in")


# ---------- Genus Abundance Bar Plot ---------- #
gen_plot <- plot_bar(phylo, x="sample", fill="genus") + 
  #facet_wrap(~Treatment, scales="free_x", nrow = 1) + 
  scale_fill_viridis_d(option="turbo") +
  labs(x = "Samples") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x =element_text(angle=90),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.grid = element_blank()
  )

#gen_plot

ggsave(gen_plot, file = file.path(exportDir, "genus_abundence.png"), height = 8, width = 8, units = "in")

# THE END #
###########
