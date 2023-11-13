library(ggplot2)
library(plyr)
library(dplyr)
library(ggpubr)


MERGED.TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/02_reduced_designs/reduced_designs_experiment_results_merged.tsv"
OUTPUT.DIR <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/figures"
DODGE.WIDTH <- 0.5
ERRORBAR.WIDTH <- 0.384
BOXPLOT.WIDTH <- 0.618
LINE.CIRCLE.SIZE <- 2.62
SOLVER.COLORS <- c("ACE" = "#D65DB1",
                   "Random" = "#FF9671", # F05B61
                   "Repeated" = "#ABABAA")
COMMON.THEME <- theme(axis.title = element_text(size = 14),
                      axis.text = element_text(size = 12),
                      legend.title = element_text(size = 12),
                      legend.text = element_text(size = 12))
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
    scale_x_continuous(breaks = seq(1,15,1)) +
    scale_y_continuous(limits = c(num.first.round.pools, ymax), breaks = breaks) +
    theme_pubr() + COMMON.THEME
  return(plot)
}
plot.200.10.3.precision <- PlotPrecision(df.plot = df.plot, group = "200/10/3")
plot.400.20.3.precision <- PlotPrecision(df.plot = df.plot, group = "400/20/3")
plot.100.10.3.num.pools <- PlotNumTotalPools(df.plot = df.plot, group = "100/10/3", ymax = 130, breaks = c(30,40,70,100,130))
plot.200.10.3.num.pools <- PlotNumTotalPools(df.plot = df.plot, group = "200/10/3", ymax = 260, breaks = c(60,80,140,200,260))
plot.400.20.3.num.pools <- PlotNumTotalPools(df.plot = df.plot, group = "400/20/3", ymax = 460, breaks = c(60,140,220,300,380,460))
plot.800.25.3.num.pools <- PlotNumTotalPools(df.plot = df.plot, group = "800/25/3", ymax = 900, breaks = c(96,300,600,900))

figure <- ggarrange(plotlist = list(
  plot.200.10.3.precision, 
  plot.400.20.3.precision,
  plot.200.10.3.num.pools,
  plot.400.20.3.num.pools,
  plot.100.10.3.num.pools,
  plot.800.25.3.num.pools
  ),
  ncol = 2, 
  nrow = 3, 
  align = "hv", 
  common.legend = TRUE, 
  heights = c(1,2,2)
)
ggsave(plot = figure, filename = paste0(OUTPUT.DIR, "/supplementary_figure_1.pdf"),
       width = 16, height = 16, dpi = 300)
