library(ggplot2)
library(ggpubr)
library(plyr)


# Step 1. Define constants
OUTPUT.DIR <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/figures"
DODGE.WIDTH <- 0.5
ERRORBAR.WIDTH <- 0.146
BOXPLOT.WIDTH <- 0.263
COMMON.THEME <- theme(axis.title = element_text(size = 14),
                      axis.text = element_text(size = 12),
                      legend.title = element_text(size = 12),
                      legend.text = element_text(size = 12))


# Step 2. Plot panel A
CSV.FILE.1 <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/seqsim_model/outgroup_euclidean_not_finetuned.csv"
CSV.FILE.2 <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/seqsim_model/ingroup_euclidean_not_finetuned.csv"
CSV.FILE.3 <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/seqsim_model/outgroup_euclidean_finetuned.csv"
CSV.FILE.4 <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/seqsim_model/ingroup_euclidean_finetuned.csv"
OUTGROUP.HEX <- "#173ccc"
INGROUP.HEX <- "#e1151b"
FILL.COLORS <- c("ingroup" = INGROUP.HEX,
                 "outgroup" = OUTGROUP.HEX)
df.outgroup.not.finetuned <- read.csv(CSV.FILE.1)
df.outgroup.not.finetuned$group <- "outgroup"
df.ingroup.not.finetuned <- read.csv(CSV.FILE.2)
df.ingroup.not.finetuned$group <- "ingroup"
df.not.finetuned <- rbind(df.outgroup.not.finetuned, df.ingroup.not.finetuned)
df.outgroup.finetuned <- read.csv(CSV.FILE.3)
df.outgroup.finetuned$group <- "outgroup"
df.ingroup.finetuned <- read.csv(CSV.FILE.4)
df.ingroup.finetuned$group <- "ingroup"
df.finetuned <- rbind(df.outgroup.finetuned, df.ingroup.finetuned)
plot.not.finedtuned <- ggplot(df.not.finetuned, aes(x = X0, fill = group)) +
  geom_density(alpha = 0.618) +
  xlab("Euclidean Similarity") + ylab("Density") +
  scale_x_continuous(limits = c(0,1)) +
  scale_fill_manual(values = FILL.COLORS) +
  theme_pubr() + COMMON.THEME
plot.finedtuned <- ggplot(df.finetuned, aes(x = X0, fill = group)) +
  geom_density(alpha = 0.618) +
  xlab("Euclidean Similarity") + ylab("Density") +
  scale_x_continuous(limits = c(0,1)) +
  scale_fill_manual(values = FILL.COLORS) +
  theme_pubr() + COMMON.THEME
figure.panel.a <- ggarrange(plotlist = list(plot.not.finedtuned, plot.finedtuned),
                            ncol = 2, align = "hv")
print(figure.panel.a)


# Step 3. Plot panel B
MERGED.TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/05_alanine_scanning/alanine_scanning_experiment_results_merged.tsv"
SOLVER.COLORS <- c("ACE" = "#845EC2",
                   "ACE (w/o clustering)" = "#D65DB1",
                   "Randomized" = "#FF9671", # F05B61
                   "Repeated" = "#ABABAA")
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
         "Randomized",
         "Repeated")
)
df.plot$solver <- factor(
  df.plot$solver,
  levels = c("ACE",
             "ACE (w/o clustering)",
             "Randomized",
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
      "1","5","10","15","20"
    )
  )
  plot <- ggplot(df.plot.temp, aes(x = perc_positive_peptide_sequences, y = precision, fill = solver)) +
    geom_col(width = 0.8, colour = "black", position = "dodge") +
    geom_errorbar(aes(ymin = precision - sd, ymax = ifelse(precision + sd > 1, 1, precision + sd)), width = 0.2, position = position_dodge(0.8)) +
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
    geom_line(size = 0.8) +
    geom_hline(yintercept = num.first.round.pools,  linetype = "dashed") +
    geom_errorbar(aes(ymin = predicted_total_pools - sd, ymax = predicted_total_pools + sd), width = 0.384, size = 0.8) +
    geom_point(size = 4) +
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


# Step 4. Plot panel C
MERGED.TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/06_real_world_datasets/tan_et_al_inf_gen_evol_2021/177peptides_10perpool_3x/tan_et_al_inf_gen_evol_2021_experiment_results_582168310.tsv"
FILL.COLORS <- c("ACE (w/ Extra Pools)" = "#845EC2",
                 "ACE (w/o Clustering, w/ Extra Pools)" = "#D65DB1",
                 "Randomized Block Design" = "#FF9671", # F05B61
                 "Repeated Block Design" = "#ABABAA")
df.plot <- read.csv(MERGED.TSV.FILE, sep = '\t')
df.plot$group <- paste0(df.plot$num_peptides, "/", df.plot$num_peptides_per_pool, "/", df.plot$num_coverage)
df.plot <- df.plot[df.plot$solver %in% c("ace_golfy_clusteron_extrapools",
                                         "ace_golfy_clusteroff_extrapools",
                                         "randomized_block_assignment",
                                         "repeated_block_assignment"),]
df.plot$solver <- mapvalues(
  x = df.plot$solver, 
  from = c("ace_golfy_clusteron_extrapools",
           "ace_golfy_clusteroff_extrapools",
           "randomized_block_assignment",
           "repeated_block_assignment"),
  to = c("ACE (w/ Extra Pools)",
         "ACE (w/o Clustering, w/ Extra Pools)",
         "Randomized Block Design",
         "Repeated Block Design")
)
df.plot$solver <- factor(
  df.plot$solver,
  levels = c("ACE (w/ Extra Pools)",
             "ACE (w/o Clustering, w/ Extra Pools)",
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
  stat_boxplot(geom = "errorbar", width = 0.4) +
  geom_boxplot(outlier.shape = NA, width = 0.8) +
  xlab("") + ylab("Number of Total Pools") + 
  scale_fill_manual(values = FILL.COLORS) +
  scale_y_continuous(limits = c(50, 200), breaks = c(54,100,150,200)) +
  theme_pubr() + COMMON.THEME
df.plot.temp <- df.plot %>%
  dplyr::group_by(solver) %>%
  dplyr::summarise(
    sd = sd(precision_empirical),
    precision = mean(precision_empirical)
  )
plot.precision <- ggplot(df.plot.temp, aes(x = solver, y = precision, fill = solver)) +
  geom_col(width = 0.8, colour = "black", position = "dodge") +
  geom_errorbar(aes(ymin = precision - sd, ymax = ifelse(precision + sd > 1, 1, precision + sd)), width = 0.2, position = position_dodge(0.5)) +
  xlab("") + ylab("Precision") +
  scale_y_continuous(limits = c(0,1), breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) +
  scale_fill_manual(values = FILL.COLORS) +
  theme_pubr() + COMMON.THEME


# Step 5. Plot panel D
MERGED.TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/06_real_world_datasets/cameron_et_al_sci_trans_med_2013/36peptides_6perpool_3x/cameron_et_al_sci_trans_med_2013_experiment_results_207124847.tsv"
FILL.COLORS <- c("ACE" = "#845EC2",
                 "ACE (w/o clustering)" = "#D65DB1",
                 "Randomized Block Design" = "#FF9671", # F05B61
                 "Repeated Block Design" = "#ABABAA")
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
  stat_boxplot(geom = "errorbar", width = 0.4) +
  geom_boxplot(outlier.shape = NA, width = 0.8) +
  xlab("") + ylab("Number of Total Pools") + 
  scale_fill_manual(values = FILL.COLORS) +
  scale_y_continuous(limits = c(18, 36), breaks = c(18,24,30,36)) +
  theme_pubr() + COMMON.THEME
df.plot.temp <- df.plot %>%
  dplyr::group_by(solver) %>%
  dplyr::summarise(
    sd = sd(precision_empirical),
    precision = mean(precision_empirical)
  )
plot.precision <- ggplot(df.plot.temp, aes(x = solver, y = precision, fill = solver)) +
  geom_col(width = 0.8, colour = "black", position = "dodge") +
  geom_errorbar(aes(ymin = precision - sd, ymax = ifelse(precision + sd > 1, 1, precision + sd)), width = 0.2, position = position_dodge(0.5)) +
  xlab("") + ylab("Precision") +
  scale_y_continuous(limits = c(0,1), breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) +
  scale_fill_manual(values = FILL.COLORS) +
  theme_pubr() + COMMON.THEME


# Step 6. Plot panel E
MERGED.TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/06_real_world_datasets/tiling_window_covid_spike_9mers/1265peptides_20perpool_3x/tiling_window_covid_spike_9mers_experiment_results_merged_20perpool.tsv"
FILL.COLORS <- c("ACE" = "#845EC2",
                 "ACE (w/o Clustering)" = "#D65DB1",
                 "Randomized Block Design" = "#FF9671", # F05B61
                 "Repeated Block Design" = "#ABABAA")
df.plot <- read.csv(MERGED.TSV.FILE, sep = '\t')
df.plot$group <- paste0(df.plot$num_peptides, "/", df.plot$num_peptides_per_pool, "/", df.plot$num_coverage)
df.plot$solver <- mapvalues(
  x = df.plot$solver, 
  from = c("ace_golfy_clusteron_noextrapools",
           "ace_golfy_clusteroff_noextrapools",
           "randomized_block_assignment",
           "repeated_block_assignment"),
  to = c("ACE",
         "ACE (w/o Clustering)",
         "Randomized Block Design",
         "Repeated Block Design")
)
df.plot$solver <- factor(
  df.plot$solver,
  levels = c("ACE",
             "ACE (w/o Clustering)",
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
  stat_boxplot(geom = "errorbar", width = 0.4) +
  geom_boxplot(outlier.shape = NA, width = 0.8) +
  xlab("") + ylab("Number of Total Pools") + 
  scale_fill_manual(values = FILL.COLORS) +
  theme_pubr() + COMMON.THEME
df.plot.temp <- df.plot %>%
  dplyr::group_by(solver) %>%
  dplyr::summarise(
    sd = sd(precision_empirical),
    precision = mean(precision_empirical)
  )
plot.precision <- ggplot(df.plot.temp, aes(x = solver, y = precision, fill = solver)) +
  geom_col(width = 0.8, colour = "black", position = "dodge") +
  geom_errorbar(aes(ymin = precision - sd, ymax = ifelse(precision + sd > 1, 1, precision + sd)), width = 0.2, position = position_dodge(0.5)) +
  xlab("") + ylab("Precision") +
  scale_y_continuous(limits = c(0,1), breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) +
  scale_fill_manual(values = FILL.COLORS) +
  theme_pubr() + COMMON.THEME

