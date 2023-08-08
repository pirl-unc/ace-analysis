library(ggplot2)
library(plyr)
library(dplyr)
library(ggpubr)
library(ggsignif)


MERGED.TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/alanine_scanning/alanine_scanning_experiment_results_merged.tsv"
OUTPUT.DIR <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/alanine_scanning/"
DODGE.WIDTH <- 0.5
ERRORBAR.WIDTH <- 0.146
BOXPLOT.WIDTH <- 0.8
SOLVER.COLORS <- c("Golfy (+Clustering / -Extra Pools)" = "#845EC2",
                   "Golfy (-Clustering / -Extra Pools)" = "#D65DB1",
                   "Randomized Block Design" = "#FF9671", # F05B61
                   "Repeated Block Design" = "#ABABAA")

# SOLVER.COLORS <- c("Golfy (+Clustering / -Extra Pools)" = "#1E6894",
#                    "Golfy (-Clustering / -Extra Pools)" = "#2B95D4",
#                    "Randomized Block Design" = "#DD63A1", # F05B61
#                    "Repeated Block Design" = "#ABABAA")

# Step 2. Load data
df.plot <- read.csv(MERGED.TSV.FILE, sep = '\t')

# Step 3. Postprocess data
df.plot$group <- paste0(df.plot$num_peptides, "/", df.plot$num_peptides_per_pool, "/", df.plot$num_coverage)
df.plot$perc_positive_peptide_sequences <- floor((df.plot$num_positive_peptide_sequences / df.plot$num_peptides) * 100)
df.plot <- df.plot[df.plot$solver %in% c("ace_golfy_clusteron_noextrapools",
                                         "ace_golfy_clusteroff_noextrapools",
                                         "randomized_block_assignment",
                                         "repeated_block_assignment"),]

# Step 4. Enforce ordering
df.plot$solver <- mapvalues(
  x = df.plot$solver, 
  from = c("ace_golfy_clusteron_noextrapools",
           "ace_golfy_clusteroff_noextrapools",
           "randomized_block_assignment",
           "repeated_block_assignment"),
  to = c("Golfy (+Clustering / -Extra Pools)",
         "Golfy (-Clustering / -Extra Pools)",
         "Randomized Block Design",
         "Repeated Block Design")
)
df.plot$solver <- factor(
  df.plot$solver,
  levels = c("Golfy (+Clustering / -Extra Pools)",
             "Golfy (-Clustering / -Extra Pools)",
             "Randomized Block Design",
             "Repeated Block Design")
)

# Step 5. Plot
for (group in unique(df.plot$group)) {
  df.plot.temp <- df.plot[df.plot$group == group,]
  num.peptides <- df.plot.temp$num_peptides[1]
  num.peptides.per.pool <- df.plot.temp$num_peptides_per_pool[1]
  num.coverage <- df.plot.temp$num_coverage[1]
  num.first.round.pools <- (df.plot.temp$num_peptides[1] / df.plot.temp$num_peptides_per_pool[1]) * df.plot.temp$num_coverage[1]
  num.max.pools <- max(df.plot.temp$predicted_total_pools)
  
  df.plot.temp.2 <- df.plot.temp %>%
    dplyr::group_by(solver, perc_positive_peptide_sequences) %>%
    dplyr::summarise(
      sd = sd(predicted_total_pools, na.rm = T),
      predicted_total_pools = mean(predicted_total_pools)
    )
  plot.num.total.pools <- ggplot(df.plot.temp.2, aes(perc_positive_peptide_sequences, predicted_total_pools, colour = solver)) +
    geom_line(size = 0.8) +
    geom_hline(yintercept = num.first.round.pools,  linetype = "dashed") +
    geom_errorbar(aes(ymin = predicted_total_pools - sd, ymax = predicted_total_pools + sd), width = 0.384, size = 0.8) +
    geom_point(size = 4) +
    xlab("Percentage of Positive Peptides") + ylab("Number of Total Pools") + ggtitle(group) +
    scale_color_manual(values = SOLVER.COLORS) +
    scale_x_continuous(limits = c(0.5,20.5), breaks = c(1,5,10,15,20)) +
    scale_y_continuous(limits = c(num.first.round.pools, num.max.pools * 1.05)) +
    theme_pubr() +
    theme(plot.title = element_text(size = 12, hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12))
  print(plot.num.total.pools)
  # Plot number of total pools
  # plot.num.total.pools <- ggplot(df.plot.temp, aes(x = solver, y = predicted_total_pools, fill = solver)) +
  #   geom_hline(yintercept = num.first.round.pools, linetype = "dashed") +
  #   stat_boxplot(geom = "errorbar", width = ERRORBAR.WIDTH, position = position_dodge(width = DODGE.WIDTH)) +
  #   geom_boxplot(position = position_dodge(width = DODGE.WIDTH), width = BOXPLOT.WIDTH, outlier.shape = NA) +
  #   facet_wrap(~ perc_positive_peptide_sequences, nrow =  1, strip.position = "bottom") +
  #   xlab("Percentage of Positive Peptides") + ylab("Number of Total Pools") + ggtitle(group) +
  #   scale_fill_manual(values = SOLVER.COLORS) +
  #   theme_bw() +
  #   theme(plot.title = element_text(size = 12, hjust = 0.5),
  #         axis.title = element_text(size = 12),
  #         axis.text.x = element_text(size = 0),
  #         axis.text.y = element_text(size = 12),
  #         axis.ticks.x = element_blank(),
  #         panel.grid.major.x = element_blank(),
  #         panel.grid.minor.y = element_blank(),
  #         legend.title = element_text(size = 12),
  #         legend.text = element_text(size = 12),
  #         strip.text = element_text(size = 12),
  #         strip.background = element_rect(fill = "white"))
  ggsave(plot = plot.num.total.pools, 
         filename = paste0(OUTPUT.DIR, '/total_pools_', num.peptides, 'peptides_', num.peptides.per.pool, 'perpool_', num.coverage, 'x.pdf'),
         width = 12, height = 8, dpi = 300
  )
  
  # Plot precision
  df.plot.temp.3 <- df.plot.temp %>%
    dplyr::group_by(solver, perc_positive_peptide_sequences) %>%
    dplyr::summarise(
      sd = sd(precision_empirical),
      precision = mean(precision_empirical)
    )
  df.plot.temp.3$perc_positive_peptide_sequences <- as.character(df.plot.temp.3$perc_positive_peptide_sequences)
  df.plot.temp.3$perc_positive_peptide_sequences <- factor(
    df.plot.temp.3$perc_positive_peptide_sequences,
    levels = c(
      "1","5","10","15","20"
    )
  )
  plot.precision <- ggplot(df.plot.temp.3, aes(x = perc_positive_peptide_sequences, y = precision, fill = solver)) +
    geom_col(width = 0.8, colour = "black", position = "dodge") +
    geom_errorbar(aes(ymin = precision - sd, ymax = ifelse(precision + sd > 1, 1, precision + sd)), width = 0.2, position = position_dodge(0.8)) +
    xlab("Percentage of Positive Peptides") + ylab("Precision") + ggtitle(group) +
    scale_y_continuous(limits = c(0,1), breaks = c(0.00, 0.25, 0.50, 0.75, 1.00), expand = c(0,0)) +
    scale_fill_manual(values = SOLVER.COLORS) +
    theme_pubr() +
    theme(plot.title = element_text(size = 12, hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          strip.text = element_text(size = 12),
          strip.background = element_rect(fill = "white"))
  print(plot.precision)
  ggsave(plot = plot.precision, 
         filename = paste0(OUTPUT.DIR, '/precision_', num.peptides, 'peptides_', num.peptides.per.pool, 'perpool_', num.coverage, 'x.pdf'),
         width = 8, height = 3, dpi = 300
  )
}








library(ggplot2)
library(ggsignif)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(plyr)
library(dplyr)
library(ggbreak)
library(ggprism)


# Step 1. Define constants
MERGED.TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/alanine_scanning/alanine_scanning_experiment_results_merged.tsv"
OUTPUT.DIR <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/alanine_scanning"
DODGE.WIDTH <- 0.5
ERRORBAR.WIDTH <- 0.2
BOXPLOT.WIDTH <- 0.4
SOLVER.COLORS <- c("Golfy (-Clustering / -Extra Pools)" = "#abd5ee",
                   "Randomized Block Design" = "#BDBDBD",
                   "Repeated Block Design" = "#d7d7d7")

# Step 2. Load data
df.plot <- read.csv(MERGED.TSV.FILE, sep = '\t')

# Step 3. Postprocess data
df.plot$group <- paste0(df.plot$num_peptides, "/", df.plot$num_peptides_per_pool, "/", df.plot$num_coverage)
df.plot$perc_positive_peptide_sequences <- floor((df.plot$num_positive_peptide_sequences / df.plot$num_peptides) * 100)
df.plot$perc_positive_peptide_sequences <- as.character(df.plot$perc_positive_peptide_sequences)

# Step 4. Enforce ordering
df.plot$perc_positive_peptide_sequences <- factor(
  df.plot$perc_positive_peptide_sequences,
  levels = c(
    "1",
    "5",
    "10",
    "15",
    "20"
  )
)
df.plot$solver <- factor(
  df.plot$solver,
  levels = c(
    "ace_golfy_clusteron_noextrapools",
    "ace_golfy_clusteroff_noextrapools",
    "randomized_block_assignment",
    "repeated_block_assignment"
    )
)
df.plot$group <- factor(
  df.plot$group,
  levels = c(
    "90/9/3",
    "180/9/3",
    "360/18/3",
    "720/18/3"
  )
)
df.plot$solver <- mapvalues(
  x = df.plot$solver, 
  from = c("ace_golfy_clusteron_noextrapools",
           "ace_golfy_clusteroff_noextrapools",
           "randomized_block_assignment",
           "repeated_block_assignment"),
  to = c("Golfy (+Clustering / -Extra Pools)",
         "Golfy (-Clustering / -Extra Pools)",
         "Randomized Block Design",
         "Repeated Block Design")
)
df.plot$group <- mapvalues(
  x = df.plot$group, 
  from = c("90/18/3",
           "180/18/3",
           "360/18/3",
           "720/18/3"),
  to = c("90 Peptides | 9 Peptides Per Pool | 3 Replicates",
         "180 Peptides | 9 Peptides Per Pool | 3 Replicates",
         "360 Peptides | 18 Peptides Per Pool | 3 Replicates",
         "720 Peptides | 18 Peptides Per Pool | 3 Replicates")
)

# Step 5. Plot
for (group in unique(df.plot$group)) {
  df.plot.temp <- df.plot[df.plot$group == group,]
  num.peptides <- df.plot.temp$num_peptides[1]
  num.peptides.per.pool <- df.plot.temp$num_peptides_per_pool[1]
  num.coverage <- df.plot.temp$num_coverage[1]
  num.first.round.pools <- (df.plot.temp$num_peptides[1] / df.plot.temp$num_peptides_per_pool[1]) * df.plot.temp$num_coverage[1]
  num.naive.approach.total.pools <- df.plot.temp$num_peptides[1] * df.plot.temp$num_coverage[1]
  max.total.pools <- max(df.plot.temp$predicted_total_pools)
  min.total.pools <- min(df.plot.temp$predicted_total_pools)

  # Step 1. Compute statistical significance
  # df.p.val <- df.plot.temp %>%
  #   dplyr::filter(perc_positive_peptide_sequences != 1) %>%
  #   rstatix::group_by(perc_positive_peptide_sequences) %>%
  #   rstatix::t_test(predicted_total_pools ~ solver) %>%
  #   rstatix::add_significance(p.col = "p") %>%
  #   rstatix::add_xy_position(x = "perc_positive_peptide_sequences", dodge = DODGE.WIDTH) # important for positioning!
  
  # Step 2. Plot number of total pools
  plot.num.total.pools <- ggplot(df.plot.temp, aes(x = factor(perc_positive_peptide_sequences),
                                                   y = predicted_total_pools)) +
    xlab("Percentage of Positive Peptides") + ylab("Number of Total Pools") + ggtitle(group) +
    geom_hline(yintercept = num.first.round.pools, linetype = "dashed") +
    stat_boxplot(mapping = aes(fill = solver), geom = "errorbar", width = ERRORBAR.WIDTH, position = position_dodge(width = DODGE.WIDTH)) +
    geom_boxplot(aes(fill = solver), position = position_dodge(width = DODGE.WIDTH), width = BOXPLOT.WIDTH, outlier.shape = NA) +
    # add_pvalue(data = df.p.val,
    #            xmin = "xmin",
    #            xmax = "xmax",
    #            label = "{p.signif}",
    #            tip.length = 0) +
    scale_fill_manual(values = SOLVER.COLORS) +
    theme_bw() +
    theme(plot.title = element_text(size = 12, hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          strip.text = element_text(size = 12))
  ggsave(plot = plot.num.total.pools, 
         filename = paste0(OUTPUT.DIR, '/total_pools_', num.peptides, 'peptides_', num.peptides.per.pool, 'perpool_', num.coverage, 'x.pdf'),
         width = 9, height = 6, dpi = 300
  )
  
  # Step 3. Plot precision
  plot.precision <- ggplot(df.plot.temp, aes(x = factor(perc_positive_peptide_sequences),
                                             y = precision_empirical)) +
    xlab("Percentage of Positive Peptides") + ylab("Precision") + ggtitle(group) +
    stat_boxplot(mapping = aes(fill = solver), geom = "errorbar", width = ERRORBAR.WIDTH, position = position_dodge(width = DODGE.WIDTH)) +
    geom_boxplot(aes(fill = solver), position = position_dodge(width = DODGE.WIDTH), width = BOXPLOT.WIDTH, outlier.shape = NA) +
    # add_pvalue(data = df.p.val,
    #            xmin = "xmin",
    #            xmax = "xmax",
    #            label = "{p.signif}",
    #            tip.length = 0) +
    scale_fill_manual(values = SOLVER.COLORS) +
    theme_bw() +
    theme(plot.title = element_text(size = 12, hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          strip.text = element_text(size = 12))
  ggsave(plot = plot.precision, 
         filename = paste0(OUTPUT.DIR, '/precision_', num.peptides, 'peptides_', num.peptides.per.pool, 'perpool_', num.coverage, 'x.pdf'),
         width = 9, height = 6, dpi = 300
  )
}


