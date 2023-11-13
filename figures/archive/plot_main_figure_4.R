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
plot.empty <- ggplot()


# Step 2. Plot panels A and B
INGROUP.CSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/model_snapshot/ingroup_euclidean_sims.csv"
OUTGROUP.CSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/model_snapshot/outgroup_euclidean_sims.csv"
df.ingroup <- read.csv(INGROUP.CSV.FILE)
df.outgroup <- read.csv(OUTGROUP.CSV.FILE)
df.not.finetuned <- rbind(
  data.frame(
    x = df.ingroup$ingroup_euclidean,
    group = rep("ingroup", nrow(df.ingroup))
  ),
  data.frame(
    x = df.outgroup$outgroup_euclidean,
    group = rep("outgroup", nrow(df.outgroup))
  )
)
df.finetuned <- rbind(
  data.frame(
    x = df.ingroup$fine_tuned_euclidean,
    group = rep("ingroup", nrow(df.ingroup))
  ),
  data.frame(
    x = df.outgroup$finetuned_outgroup_euclidean,
    group = rep("outgroup", nrow(df.outgroup))
  )
)
df.levenshtein <- rbind(
  data.frame(
    x = df.ingroup$levenshtein_euclidean,
    group = rep("ingroup", nrow(df.ingroup))
  ),
  data.frame(
    x = df.outgroup$levenshtein_euclidean,
    group = rep("outgroup", nrow(df.outgroup))
  )
)
plot.not.finedtuned <- ggplot(df.not.finetuned, aes(x = x, fill = group)) +
  geom_density(alpha = 0.618) +
  xlab("Euclidean Similarity") + ylab("Density") + ggtitle("Not Fine-tuned") +
  scale_x_continuous(limits = c(0,1)) +
  scale_fill_manual(values = GROUP.COLORS) +
  theme_pubr() + COMMON.THEME
plot.finedtuned <- ggplot(df.finetuned, aes(x = x, fill = group)) +
  geom_density(alpha = 0.618) +
  xlab("Euclidean Similarity") + ylab("Density") + ggtitle("Fine-tuned") +
  scale_x_continuous(limits = c(0,1)) +
  scale_fill_manual(values = GROUP.COLORS) +
  theme_pubr() + COMMON.THEME
plot.levenshtein <- ggplot(df.levenshtein, aes(x = x, fill = group)) +
  geom_density(alpha = 0.618) +
  xlab("Euclidean Similarity") + ylab("Density") + ggtitle("Levenshtein") +
  scale_x_continuous(limits = c(0,1)) +
  scale_fill_manual(values = GROUP.COLORS) +
  theme_pubr() + COMMON.THEME
figure.panel.a <- ggarrange(plotlist = list(plot.not.finedtuned,
                                            plot.finedtuned,
                                            plot.levenshtein),
                            ncol = 3, align = "hv")
print(figure.panel.a)


# Step 3. Plot panel C and D
MERGED.TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/05_alanine_scanning/held_out_data/alanine_scanning_experiment_results_merged.tsv"
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
panel.b.1 <- PlotPrecision(df.plot = df.plot, group = "90/9/3")
panel.b.2 <- PlotPrecision(df.plot = df.plot, group = "720/9/3")
panel.b.3 <- PlotNumTotalPools(df.plot = df.plot, group = "90/9/3", ymax = 120, breaks = c(30,60,90,120))
panel.b.4 <- PlotNumTotalPools(df.plot = df.plot, group = "720/9/3", ymax = 1000, breaks = c(240, 500, 750, 1000))
figure.panel.b.temp.1 <- ggarrange(plotlist = list(panel.b.1,
                                                   panel.b.2),
                                  ncol = 1, nrow = 2, align = "hv",
                                  common.legend = T)
figure.panel.b.temp.2 <- ggarrange(plotlist = list(panel.b.3,
                                                   panel.b.4),
                                   ncol = 2, nrow = 1, align = "hv",
                                   common.legend = T)
figure.panel.b <- ggarrange(plotlist = list(plot.empty,
                                            figure.panel.b.temp.1,
                                            figure.panel.b.temp.2),
                            ncol = 3, nrow = 1, align = "hv", widths = c(1,2,4))
print(figure.panel.b)


# Step 4. Plot panel E
MERGED.TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/06_real_world_datasets/tan_et_al_inf_gen_evol_2021/noextrapools/177peptides_10perpool_3x/tan_et_al_inf_gen_evol_2021_experiment_results_585570736.tsv"
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
  scale_y_continuous(limits = c(50, 200), breaks = c(54,100,150,200)) +
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
figure.panel.c <- ggarrange(plotlist = list(plot.precision, plot.num.total.pools),
                            ncol = 2, align = "hv", common.legend = T)
print(figure.panel.c)


# Step 5. Plot panel F
MERGED.TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/06_real_world_datasets/cameron_et_al_sci_trans_med_2013/36peptides_6perpool_3x/cameron_et_al_sci_trans_med_2013_experiment_results_207124847.tsv"
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
  scale_y_continuous(limits = c(18, 36), breaks = c(18,24,30,36)) +
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
figure.panel.d <- ggarrange(plotlist = list(plot.precision, plot.num.total.pools),
                            ncol = 2, align = "hv", common.legend = T)
print(figure.panel.d)
figure.panels.c.d <- ggarrange(plotlist = list(
  figure.panel.c, 
  figure.panel.d,
  plot.empty
), ncol = 3, align = "hv", common.legend = T, widths = c(4,4,3))
print(figure.panels.c.d)


# Step 6. Plot figure 6.
figure.6 <- ggarrange(plotlist = list(
  figure.panel.a,
  figure.panel.b,
  figure.panels.c.d
), ncol = 1, align = "hv", heights = c(2,4,3))
print(figure.6)
ggsave(plot = figure.6, filename = paste0(OUTPUT.DIR, "/main_figure_4.pdf"),
       width = 16, height = 16, dpi = 300)


