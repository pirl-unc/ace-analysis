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


# Step 2. Plot panel A
CSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/01_precision_recall_violations_iterations/precision_recall_violations.csv"
df <- read.csv(CSV.FILE)
df.plot <- df %>%
  dplyr::group_by(iteration) %>%
  dplyr::summarise(precision.mean = mean(precision),
            precision.std = sd(precision),
            recall.mean = mean(recall),
            recall.std = sd(recall),
            f1.mean = mean(f1),
            f1.std = sd(f1),
            violations.mean = mean(num_violations),
            violations.std = sd(num_violations))
df.plot.precision <- df.plot %>%
  dplyr::reframe(iteration = iteration,
                 precision.mean = precision.mean,
                 precision.lower = ifelse(precision.mean - precision.std < 0, 0, precision.mean - precision.std),
                 precision.upper = ifelse(precision.mean + precision.std > 1, 1, precision.mean + precision.std))
plot.precision <- ggplot(df.plot.precision, aes(x = iteration)) +
  geom_ribbon(aes(ymax = precision.upper, ymin = precision.lower), fill = "#fcf3f8") +
  geom_line(aes(y = precision.mean), color = "#a22867", linewidth = 0.618) +
  geom_line(aes(y = precision.lower), color = "#d1458e", linewidth = 0.618) +
  geom_line(aes(y = precision.upper), color = "#e08284", linewidth = 0.618) +
  xlab("Iteration") + ylab("Precision") + ggtitle("Precision") +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1.0)) +
  theme_pubr() + COMMON.THEME
df.plot.recall <- df.plot %>%
  dplyr::reframe(iteration = iteration,
                 recall.mean = recall.mean,
                 recall.lower = ifelse(recall.mean - recall.std < 0, 0, recall.mean - recall.std),
                 recall.upper = ifelse(recall.mean + recall.std > 1, 1, recall.mean + recall.std))
plot.recall <- ggplot(df.plot.recall, aes(x = iteration)) +
  geom_ribbon(aes(ymax = recall.upper, ymin = recall.lower), fill = "#fcf3f8") +
  geom_line(aes(y = recall.mean), color = "#a22867", linewidth = 0.618) +
  geom_line(aes(y = recall.lower), color = "#d1458e", linewidth = 0.618) +
  geom_line(aes(y = recall.upper), color = "#e08284", linewidth = 0.618) +
  xlab("Iteration") + ylab("Sensitivity") + ggtitle("Sensitivity") +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1.0)) +
  theme_pubr() + COMMON.THEME
df.plot.violations <- df.plot %>%
  dplyr::reframe(iteration = iteration,
                 violations.mean = violations.mean,
                 violations.lower = ifelse(violations.mean - violations.std < 0, 0, violations.mean - violations.std),
                 violations.upper = ifelse(violations.mean + violations.std > 1, 1, violations.mean + violations.std))
plot.violations <- ggplot(df.plot.violations, aes(x = iteration)) +
  geom_line(aes(y = violations.mean), color = "#a22867", linewidth = 0.618) +
  xlab("Iteration") + ylab("Number of Violations") + ggtitle("Violations") +
  scale_y_continuous(limits = c(0,8000), breaks = c(0,2000,4000,6000,8000)) +
  theme_pubr() + COMMON.THEME
figure.panel.a <- ggarrange(plotlist = list(plot.precision, plot.recall, plot.violations),
                            ncol = 1, align = "hv")
print(figure.panel.a)


# Step 2. Plot panel B
MERGED.TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/02_reduced_designs/reduced_designs_experiment_results_merged.tsv"
SOLVER.COLORS <- c("ACE" = "#D65DB1",
                   "Random" = "#FF9671", # F05B61
                   "Repeated" = "#ABABAA")
df.plot <- read.csv(MERGED.TSV.FILE, sep = '\t')
df.plot$group <- paste0(df.plot$num_peptides, "/", df.plot$num_peptides_per_pool, "/", df.plot$num_coverage)
df.plot$perc_positive_peptide_sequences <- floor((df.plot$num_positive_peptide_sequences / df.plot$num_peptides) * 100)
df.plot <- df.plot[df.plot$solver %in% c("ace_golfy_clusteroff_noextrapools",
                                         "randomized_block_assignment",
                                         "repeated_block_assignment"),]
df.plot$solver <- mapvalues(
  x = df.plot$solver, 
  from = c("ace_golfy_clusteroff_noextrapools",
           "randomized_block_assignment",
           "repeated_block_assignment"),
  to = c("ACE",
         "Random",
         "Repeated")
)
df.plot$solver <- factor(
  df.plot$solver,
  levels = c("ACE",
             "Random",
             "Repeated")
)
PlotPrecision <- function(df.plot, group) {
  df.plot.temp <- df.plot[df.plot$group == group,]
  df.plot.temp <- df.plot.temp %>%
    dplyr::group_by(solver, perc_positive_peptide_sequences) %>%
    dplyr::summarise(
      sd = sd(precision_empirical),
      precision = mean(precision_empirical)
    )
  df.plot.temp$perc_positive_peptide_sequences <- as.character(df.plot.temp$perc_positive_peptide_sequences)
  df.plot.temp$perc_positive_peptide_sequences <- factor(
    df.plot.temp$perc_positive_peptide_sequences,
    levels = c(
      "1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"
    )
  )
  plot.precision <- ggplot(df.plot.temp, aes(x = perc_positive_peptide_sequences, y = precision, fill = solver)) +
    geom_col(width = BOXPLOT.WIDTH, colour = "black", position = "dodge") +
    geom_errorbar(aes(ymin = precision - sd, ymax = ifelse(precision + sd > 1, 1, precision + sd)), width = ERRORBAR.WIDTH, position = position_dodge(BOXPLOT.WIDTH)) +
    xlab("Percentage of Positive Peptides") + ylab("Precision") + ggtitle(group) +
    scale_y_continuous(limits = c(0,1), breaks = c(0.00, 0.25, 0.50, 0.75, 1.00), expand = c(0,0)) +
    scale_fill_manual(values = SOLVER.COLORS) +
    theme_pubr() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12)) +
    COMMON.THEME
  return(plot.precision)
}
plot.100.10.3 <- PlotPrecision(df.plot = df.plot, group = "100/10/3")
plot.800.25.3 <- PlotPrecision(df.plot = df.plot, group = "800/25/3")
figure.panel.b <- ggarrange(plotlist = list(plot.100.10.3, plot.800.25.3),
                            ncol = 1, align = "hv", common.legend = TRUE, legend = "right")
print(figure.panel.b)


# Step 3. Plot panel C
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
plot.100.10.3 <- PlotAUCROC(df.plot = df.plot, group = "100/10/3")
plot.800.25.3 <- PlotAUCROC(df.plot = df.plot, group = "800/25/3")
figure.panel.c <- ggarrange(plotlist = list(plot.100.10.3, plot.800.25.3),
                            ncol = 2, align = "hv", common.legend = TRUE, legend = "bottom")
print(figure.panel.c)


# Step 4. Generate figure 1
figure.1.top <- ggarrange(plotlist = list(figure.panel.a, figure.panel.b), 
                          ncol = 2, widths = c(1,2), align = "hv", legend = "bottom")
figure.1 <- ggarrange(plotlist = list(figure.1.top, figure.panel.c), common.legend = TRUE,
                      ncol = 1, align = "hv", heights = c(1,1), legend = "bottom")
print(figure.1)
ggsave(plot = figure.1, filename = paste0(OUTPUT.DIR, "/main_figure_3.pdf"),
       width = 16, height = 12, dpi = 300)
