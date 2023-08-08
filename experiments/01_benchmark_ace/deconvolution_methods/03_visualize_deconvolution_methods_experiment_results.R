library(ggplot2)
library(dplyr)
library(ggpubr)


# Step 1. Define constants
MERGED.TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/deconvolution_methods/deconvolution_methods_experiment_results_merged.tsv"
OUTPUT.DIR <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/deconvolution_methods"
HIT.POOLS.HEX <- "#8da0cb"
EM.HEX <- "#66c2a5"
LASSO.HEX <- "#fc8d62"
FILL.COLORS <- c("empirical" = HIT.POOLS.HEX,
                 "em" = EM.HEX,
                 "lasso" = LASSO.HEX)
THEME <- theme(axis.title = element_text(size = 16),
               axis.text = element_text(size = 14),
                legend.title = element_text(size = 12),
                   legend.text = element_text(size = 12))

df <- read.csv(MERGED.TSV.FILE, sep = '\t')
df$group <- paste0(df$num_peptides, "/", df$num_peptides_per_pool, "/", df$num_coverage)
df$perc_positive_peptide_sequences <- floor((df$num_positive_peptide_sequences / df$num_peptides) * 100)
df$perc_positive_peptide_sequences <- as.character(df$perc_positive_peptide_sequences)
df <- df[df$solver %in% c("ace_golfy_clusteroff_noextrapools",
                          "randomized_block_assignment",
                          "repeated_block_assignment"),]

# Step 2. Enforce ordering
df$perc_positive_peptide_sequences <- factor(
  df$perc_positive_peptide_sequences,
  levels = c(
    "1",
    "2",
    "3",
    "4",
    "5",
    "6",
    "7",
    "8",
    "9",
    "10",
    "11",
    "12",
    "13",
    "14",
    "15"
  )
)

# Step 3. Melt dataframe
HelperFn <- function(x) {
  groups <- c()
  solvers <- c()
  deconvolution.methods <- c()
  df.temp <- data.frame(
    group = rep(x[['group']], 3),
    solver = rep(x[['solver']], 3),
    deconvolution.methods = c('empirical','em','lasso'),
    perc_positive_peptide_sequences = rep(x[['perc_positive_peptide_sequences']], 3),
    auroc = c(as.numeric(x[['aucroc_score_empirical']]), 
              as.numeric(x[['aucroc_score_em']]), 
              as.numeric(x[['aucroc_score_lasso']])),
    method = c(paste0(x[['solver']], "-empirical"), paste0(x[['solver']], "-em"), paste0(x[['solver']], "-lasso"))
  )
  return(df.temp)
}
df.plot <- apply(df, MARGIN = 1, HelperFn)
df.plot <- bind_rows(df.plot)
df.plot <- df.plot %>%
  dplyr::filter(solver == 'ace_golfy_clusteroff_noextrapools')

# Step 4. Plot
group <- "100/10/3"
for (group in unique(df.plot$group)) {
  num.peptides <- strsplit(group, split = "\\/")[[1]][1]
  num.peptides.per.pool <- strsplit(group, split = "\\/")[[1]][2]
  num.coverage <- strsplit(group, split = "\\/")[[1]][3]
  df.plot.temp <- df.plot[df.plot$group == group,]
  df.plot.temp$perc_positive_peptide_sequences <- factor(
    df.plot.temp$perc_positive_peptide_sequences,
    levels = c(
      "1","2","3","4","5",
      "6","7","8","9","10",
      "11","12","13","14","15"
    )
  )
  df.plot.temp$deconvolution.methods <- factor(
    x = df.plot.temp$deconvolution.methods,
    levels = c("em", "lasso", "empirical")
  )
  plot.auroc.score <- ggplot(df.plot.temp, aes(x = perc_positive_peptide_sequences, y = auroc, fill = deconvolution.methods)) +
    stat_boxplot(geom = "errorbar", position = position_dodge(0.7), width = 0.384) + 
    geom_boxplot(position = position_dodge(0.7), outlier.shape = NA, width = 0.618) +
    xlab("Percentage of Positive Peptides") + ylab("AUCROC") +
    scale_fill_manual(values = FILL.COLORS) +
    scale_y_continuous(limits = c(0.5,1), breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
    theme_pubr() +
    theme(plot.title = element_text(size = 12, hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12))
  ggsave(
    plot = plot.auroc.score,
    filename = paste0(
      OUTPUT.DIR, "/", num.peptides, "peptides_", num.peptides.per.pool, "perpool_", num.coverage, "x.pdf"
    ),
    width = 16, height = 9, dpi = 300
  )
}


