library(ggplot2)
library(plyr)
library(dplyr)
library(ggpubr)


# Step 1. Define constants
OUTPUT.DIR <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/figures"
ERRORBAR.WIDTH <- 0.384
BOXPLOT.WIDTH <- 0.618
COMMON.THEME <- theme(axis.title = element_text(size = 14),
                      axis.text = element_text(size = 12),
                      legend.title = element_text(size = 12),
                      legend.text = element_text(size = 12))


# Step 2. Plot panel C
MERGED.TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/03_deconvolution_methods/deconvolution_methods_experiment_results_merged.tsv"
HIT.POOLS.HEX <- "#8da0cb"
EM.HEX <- "#66c2a5"
LASSO.HEX <- "#fc8d62"
FILL.COLORS <- c("empirical" = HIT.POOLS.HEX,
                 "em" = EM.HEX,
                 "lasso" = LASSO.HEX)
df <- read.csv(MERGED.TSV.FILE, sep = '\t')
df$group <- paste0(df$num_peptides, "/", df$num_peptides_per_pool, "/", df$num_coverage)
df$perc_positive_peptide_sequences <- floor((df$num_positive_peptide_sequences / df$num_peptides) * 100)
df$perc_positive_peptide_sequences <- as.character(df$perc_positive_peptide_sequences)
df <- df[df$solver %in% c("ace_golfy_clusteroff_noextrapools",
                          "randomized_block_assignment",
                          "repeated_block_assignment"),]
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
PlotAUCROC <- function(df.plot, group) {
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
    stat_boxplot(geom = "errorbar", position = position_dodge(BOXPLOT.WIDTH), width = ERRORBAR.WIDTH) + 
    geom_boxplot(position = position_dodge(BOXPLOT.WIDTH), outlier.shape = NA, width = BOXPLOT.WIDTH) +
    xlab("Percentage of Positive Peptides") + ylab("AUCROC") + ggtitle(group) + 
    scale_fill_manual(values = FILL.COLORS) +
    scale_y_continuous(limits = c(0.5,1), breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
    theme_pubr() + COMMON.THEME
  return(plot.auroc.score)
}
plot.200.10.3 <- PlotAUCROC(df.plot = df.plot, group = "200/10/3")
plot.400.20.3 <- PlotAUCROC(df.plot = df.plot, group = "400/20/3")
figure <- ggarrange(plotlist = list(plot.200.10.3, plot.400.20.3),
                    ncol = 2, align = "hv", common.legend = TRUE, legend = "bottom")
print(figure)
ggsave(plot = figure, filename = paste0(OUTPUT.DIR, "/supplementary_figure_3.pdf"),
       width = 16, height = 8, dpi = 300)
