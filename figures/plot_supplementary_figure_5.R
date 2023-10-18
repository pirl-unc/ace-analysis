library(ggplot2)
library(ggpubr)
library(plyr)


MERGED.TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/05_alanine_scanning/training_data/alanine_scanning_experiment_results_merged.tsv"
OUTPUT.DIR <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/figures"
DODGE.WIDTH <- 0.5
ERRORBAR.WIDTH <- 0.384
BOXPLOT.WIDTH <- 0.618
LINE.CIRCLE.SIZE <- 2.62
SOLVER.COLORS <- c("ACE" = "#845EC2",
                   "ACE (w/o clustering)" = "#D65DB1",
                   "Randomized Block Design" = "#FF9671", # F05B61
                   "Repeated Block Design" = "#ABABAA")
GROUP.COLORS <- c("ingroup" = "#e1151b",
                  "outgroup" = "#173ccc")
COMMON.THEME <- theme(axis.title = element_text(size = 14),
                      axis.text = element_text(size = 12),
                      legend.title = element_text(size = 12),
                      legend.text = element_text(size = 12))
df.plot <- read.csv(MERGED.TSV.FILE, sep = '\t')
df.plot$group <- paste0(df.plot$num_peptides, "/", df.plot$num_peptides_per_pool, "/", df.plot$num_coverage)
df.plot$perc_positive_peptide_sequences <- floor((df.plot$num_positive_peptide_sequences / df.plot$num_peptides) * 100)
df.plot$solver <- mapvalues(
  x = df.plot$solver, 
  from = c("ace_golfy_clusteron_noextrapools",
           "ace_golfy_clusteroff_noextrapools",
           "randomized_block_assignment",
           "repeated_block_assignment"),
  to = c("ACE",
         "ACE (w/o clustering)",
         "Randomized Block Design",
         "Repeated Block Design")
)
df.plot$solver <- factor(
  df.plot$solver,
  levels = c("ACE",
             "ACE (w/o clustering)",
             "Randomized Block Design",
             "Repeated Block Design")
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
      "1","5","10","15","20"
    )
  )
  plot <- ggplot(df.plot.temp, aes(x = perc_positive_peptide_sequences, y = precision, fill = solver)) +
    geom_col(width = BOXPLOT.WIDTH, colour = "black", position = "dodge") +
    geom_errorbar(aes(ymin = precision - sd, ymax = ifelse(precision + sd > 1, 1, precision + sd)), width = ERRORBAR.WIDTH, position = position_dodge(BOXPLOT.WIDTH)) +
    xlab("Percentage of Positive Peptides") + ylab("Precision") + ggtitle(group) +
    scale_y_continuous(limits = c(0,1), breaks = c(0.00, 0.25, 0.50, 0.75, 1.00), expand = c(0,0)) +
    scale_fill_manual(values = SOLVER.COLORS) +
    theme_pubr() + COMMON.THEME
  return(plot)
}
PlotNumTotalPools <- function(df.plot, group, ymax, breaks) {
  df.plot.temp <- df.plot[df.plot$group == group,]
  num.peptides <- df.plot.temp$num_peptides[1]
  num.peptides.per.pool <- df.plot.temp$num_peptides_per_pool[1]
  num.coverage <- df.plot.temp$num_coverage[1]
  num.first.round.pools <- (df.plot.temp$num_peptides[1] / df.plot.temp$num_peptides_per_pool[1]) * df.plot.temp$num_coverage[1]
  
  df.plot.temp <- df.plot.temp %>%
    dplyr::group_by(solver, perc_positive_peptide_sequences) %>%
    dplyr::summarise(
      sd = sd(predicted_total_pools, na.rm = T),
      predicted_total_pools = mean(predicted_total_pools)
    )
  plot <- ggplot(df.plot.temp, aes(perc_positive_peptide_sequences, predicted_total_pools, colour = solver)) +
    geom_line(size = BOXPLOT.WIDTH) +
    geom_hline(yintercept = num.first.round.pools,  linetype = "dashed") +
    geom_errorbar(aes(ymin = predicted_total_pools - sd, ymax = predicted_total_pools + sd), width = ERRORBAR.WIDTH, size = BOXPLOT.WIDTH) +
    geom_point(size = LINE.CIRCLE.SIZE) +
    xlab("Percentage of Positive Peptides") + ylab("Number of Total Pools") + ggtitle(group) +
    scale_color_manual(values = SOLVER.COLORS) +
    scale_x_continuous(limits = c(0.5,20.5), breaks = c(1,5,10,15,20)) +
    scale_y_continuous(limits = c(num.first.round.pools, ymax), breaks = breaks) +
    theme_pubr() + COMMON.THEME
  return(plot)
}
plot.1 <- PlotPrecision(df.plot = df.plot, group = "90/9/3")
plot.2 <- PlotPrecision(df.plot = df.plot, group = "180/9/3")
plot.3 <- PlotPrecision(df.plot = df.plot, group = "360/9/3")
plot.4 <- PlotPrecision(df.plot = df.plot, group = "720/9/3")
plot.5 <- PlotNumTotalPools(df.plot = df.plot, group = "90/9/3", ymax = 120, breaks = c(30,60,90,120))
plot.6 <- PlotNumTotalPools(df.plot = df.plot, group = "180/9/3", ymax = 240, breaks = c(60,120,180,240))
plot.7 <- PlotNumTotalPools(df.plot = df.plot, group = "360/9/3", ymax = 480, breaks = c(120,240,360,480))
plot.8 <- PlotNumTotalPools(df.plot = df.plot, group = "720/9/3", ymax = 1000, breaks = c(240,480,720,960))

figure <- ggarrange(plotlist = list(
  plot.1, plot.2,
  plot.5, plot.6,
  plot.3, plot.4,
  plot.7,plot.8
), ncol = 2, nrow = 4, heights = c(1,2,1,2), align = "hv", common.legend = T, legend = "bottom")
print(figure)
ggsave(plot = figure, filename = paste0(OUTPUT.DIR, "/supplementary_figure_5.pdf"),
       width = 16, height = 16, dpi = 300)
