# The purpose of this Rscript is to visualize the results of the hit count vs candidate hits experiments.

library(ggplot2)
library(dplyr)
library(plyr)
library(reshape2)


# 01. Define constants
DATA.DIR <- '/Users/leework/Documents/Research/projects/project_ace/data/processed/02_in_silico_validation_study'
OUTPUT.DIR <- '/Users/leework/Documents/Research/projects/project_ace/data/processed/02_in_silico_validation_study'

# 02. Plot function
plotResults <- function(tsv.file) {
  tsv.file.basename <- sub('\\.tsv$', '', basename(tsv.file))
  df.results <- read.csv(tsv.file, sep = '\t')
  p <- ggplot(df.results, 
              aes(x = ground_truth_hit_peptide_ids_count, 
                  y = candidate_hit_peptide_ids_count, 
                  group = ground_truth_hit_peptide_ids_count)) +
    stat_boxplot(geom = "errorbar", width = 0.25) + 
    geom_boxplot(width = 0.5, outlier.alpha = 0.0) +
    xlab("Number of Ground Truth Peptide Hits") + ylab("Number of Candidate Peptide Hits (Additional Pools)") +
    scale_x_continuous(breaks = seq(1,max(df.results$ground_truth_hit_peptide_ids_count),1)) +
    scale_y_continuous(breaks = seq(0,200,10)) +
    coord_flip() +
    theme_bw() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16))
  ggsave(filename = paste0(OUTPUT.DIR, '/', tsv.file.basename, '.pdf'),
         plot = p, width = 16, height = 9, dpi = 300)  
}

# 03. Run plot for each configuration
for (tsv.file in Sys.glob(paste0(DATA.DIR, "/*tsv"))) {
  plotResults(tsv.file = tsv.file)
}


