library(ggplot2)
library(ggpubr)
library(plyr)


# Step 1. Define constants
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

# Step 2. Plot
MERGED.TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/06_real_world_datasets/tiling_window_covid_spike_9mers/1265peptides_20perpool_3x/tiling_window_covid_spike_9mers_experiment_results_merged_20perpool.tsv"
df.plot <- read.csv(MERGED.TSV.FILE, sep = '\t')
df.plot$group <- paste0(df.plot$num_peptides, "/", df.plot$num_peptides_per_pool, "/", df.plot$num_coverage)
df.plot <- df.plot[df.plot$solver %in% c("ace_golfy_clusteron_noextrapools",
                                         "ace_golfy_clusteroff_noextrapools",
                                         "randomized_block_assignment",
                                         "repeated_block_assignment"),]
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
num.peptides <- df.plot$num_peptides[1]
num.peptides.per.pool <- df.plot$num_peptides_per_pool[1]
num.coverage <- df.plot$num_coverage[1]
num.first.round.pools <- ceiling(df.plot$num_peptides[1] / df.plot$num_peptides_per_pool[1]) * df.plot$num_coverage[1]
num.max.pools <- max(df.plot$predicted_total_pools)
plot.num.total.pools <- ggplot(df.plot, aes(x = solver, y = predicted_total_pools, fill = solver)) + 
  geom_hline(yintercept = num.first.round.pools, linetype = "dashed") +
  stat_boxplot(geom = "errorbar", width = ERRORBAR.WIDTH) +
  geom_boxplot(outlier.shape = NA, width = BOXPLOT.WIDTH) +
  xlab("") + ylab("Number of Total Pools") + 
  scale_fill_manual(values = SOLVER.COLORS) +
  scale_y_continuous(limits = c(100, 1200), breaks = c(100,300,600,900,1200)) +
  theme_pubr() + COMMON.THEME +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1))
df.plot.temp <- df.plot %>%
  dplyr::group_by(solver) %>%
  dplyr::summarise(
    sd = sd(precision_empirical),
    precision = mean(precision_empirical)
  )
plot.precision <- ggplot(df.plot.temp, aes(x = solver, y = precision, fill = solver)) +
  geom_col(width = BOXPLOT.WIDTH, colour = "black", position = "dodge") +
  geom_errorbar(aes(ymin = precision - sd, ymax = ifelse(precision + sd > 1, 1, precision + sd)), width = ERRORBAR.WIDTH, position = position_dodge(BOXPLOT.WIDTH)) +
  xlab("") + ylab("Precision") +
  scale_y_continuous(limits = c(0,1), breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) +
  scale_fill_manual(values = SOLVER.COLORS) +
  theme_pubr() + COMMON.THEME +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1))
figure <- ggarrange(plotlist = list(plot.precision, plot.num.total.pools),
                            ncol = 2, align = "hv", common.legend = T)
print(figure)
ggsave(plot = figure, filename = paste0(OUTPUT.DIR, "/supplementary_figure_6.pdf"),
       width = 8, height = 8, dpi = 300)
