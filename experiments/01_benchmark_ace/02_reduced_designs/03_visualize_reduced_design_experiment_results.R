library(ggplot2)
library(plyr)
library(dplyr)
library(ggpubr)
library(ggsignif)


MERGED.TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/reduced_designs/reduced_designs_experiment_results_merged.tsv"
OUTPUT.DIR <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/reduced_designs/"
DODGE.WIDTH <- 0.5
ERRORBAR.WIDTH <- 0.263
BOXPLOT.WIDTH <- 0.384
SOLVER.COLORS <- c("Golfy (-Clustering / -Extra Pools)" = "#D65DB1",
                   "Randomized Block Design" = "#FF9671", # F05B61
                   "Repeated Block Design" = "#ABABAA")

# Step 2. Load data
df.plot <- read.csv(MERGED.TSV.FILE, sep = '\t')

# Step 3. Postprocess data
df.plot$group <- paste0(df.plot$num_peptides, "/", df.plot$num_peptides_per_pool, "/", df.plot$num_coverage)
df.plot$perc_positive_peptide_sequences <- floor((df.plot$num_positive_peptide_sequences / df.plot$num_peptides) * 100)
df.plot <- df.plot[df.plot$solver %in% c("ace_golfy_clusteroff_noextrapools",
                                         "randomized_block_assignment",
                                         "repeated_block_assignment"),]

# Step 4. Enforce ordering
df.plot$solver <- mapvalues(
  x = df.plot$solver, 
  from = c("ace_golfy_clusteroff_noextrapools",
           "randomized_block_assignment",
           "repeated_block_assignment"),
  to = c("Golfy (-Clustering / -Extra Pools)",
         "Randomized Block Design",
         "Repeated Block Design")
)
df.plot$solver <- factor(
  df.plot$solver,
  levels = c("Golfy (-Clustering / -Extra Pools)",
             "Randomized Block Design",
             "Repeated Block Design")
)

# Step 5. Plot
group <- "100/10/3"
for (group in unique(df.plot$group)) {
  df.plot.temp <- df.plot[df.plot$group == group,]
  num.peptides <- df.plot.temp$num_peptides[1]
  num.peptides.per.pool <- df.plot.temp$num_peptides_per_pool[1]
  num.coverage <- df.plot.temp$num_coverage[1]
  num.first.round.pools <- (df.plot.temp$num_peptides[1] / df.plot.temp$num_peptides_per_pool[1]) * df.plot.temp$num_coverage[1]

  df.plot.temp.2 <- df.plot.temp %>%
    dplyr::group_by(solver, perc_positive_peptide_sequences) %>%
    dplyr::summarise(
      sd = sd(predicted_total_pools, na.rm = T),
      predicted_total_pools = mean(predicted_total_pools)
    )
  plot.num.total.pools <- ggplot(df.plot.temp.2, aes(perc_positive_peptide_sequences, predicted_total_pools, colour = solver)) +
    geom_line(size = 1) +
    geom_errorbar(aes(ymin = predicted_total_pools - sd, ymax = predicted_total_pools + sd), width = 0.2) +
    geom_point(size = 1.618) +
    xlab("Percentage of Positive Peptides") + ylab("Number of Total Pools") + ggtitle(group) +
    scale_color_manual(values = SOLVER.COLORS) +
    scale_x_continuous(limits = c(0.5,15.5), breaks = seq(1,15,1)) +
    theme_bw() +
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
         width = 16, height = 8, dpi = 300
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
      "1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"
    )
  )
  plot.precision <- ggplot(df.plot.temp.3, aes(x = perc_positive_peptide_sequences, y = precision, fill = solver)) +
    geom_col(width = 0.5, colour = "black", position = "dodge") +
    geom_errorbar(aes(ymin = precision - sd, ymax = ifelse(precision + sd > 1, 1, precision + sd)), width = 0.2, position = position_dodge(0.5)) +
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
         width = 16, height = 4, dpi = 300
  )
}

