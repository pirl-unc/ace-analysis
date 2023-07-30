library(ggplot2)
library(plyr)


MERGED.TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/reduced_designs/reduced_designs_experiment_results_merged.tsv"
OUTPUT.DIR <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/reduced_designs/"
DODGE.WIDTH <- 0.5
ERRORBAR.WIDTH <- 0.2
BOXPLOT.WIDTH <- 0.4
SOLVER.COLORS <- c("CP-SAT Solver" = "#008000",
                   "Golfy (-Clustering / -Extra Pools)" = "#F1373D",
                   "Golfy (-Clustering / +Extra Pools)" = "#003C80",
                   "Randomized Block Design" = "#BDBDBD",
                   "Repeated Block Design" = "#d7d7d7")

# Step 2. Load data
df.plot <- read.csv(MERGED.TSV.FILE, sep = '\t')

# Step 3. Postprocess data
df.plot$group <- paste0(df.plot$num_peptides, "/", df.plot$num_peptides_per_pool, "/", df.plot$num_coverage)
df.plot$perc_positive_peptide_sequences <- floor((df.plot$num_positive_peptide_sequences / df.plot$num_peptides) * 100)
df.plot$first_round_assays <- df.plot$num_pools
df.plot$second_round_assays <- df.plot$predicted_total_pools - df.plot$first_round_assays
df.plot <- df.plot[df.plot$perc_positive_peptide_sequences %in% c(1, 5, 10, 15, 20),]
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

# df.plot$group <- factor(
#   df.plot$group,
#   levels = c(
#     "90/9/3",
#     "180/9/3",
#     "360/9/3",
#     "720/9/3"
#   )
# )
df.plot$solver <- mapvalues(
  x = df.plot$solver, 
  from = c("ace_golfy_clusteroff_noextrapools",
           "ace_golfy_clusteroff_extrapools",
           "ace_sat_solver_clusteroff",
           "randomized_block_assignment",
           "repeated_block_assignment"),
  to = c("CP-SAT Solver",
         "Golfy (-Clustering / -Extra Pools)",
         "Golfy (-Clustering / +Extra Pools)",
         "Randomized Block Design",
         "Repeated Block Design")
)
df.plot$solver <- factor(
  df.plot$solver,
  levels = c("CP-SAT Solver",
        "Golfy (-Clustering / -Extra Pools)",
        "Golfy (-Clustering / +Extra Pools)",
        "Randomized Block Design",
        "Repeated Block Design"
  )
)
# df.plot$group <- mapvalues(
#   x = df.plot$group, 
#   from = c("90/9/3",
#            "180/9/3",
#            "360/9/3",
#            "720/9/3"),
#   to = c("90 Peptides | 9 Peptides Per Pool | 3 Replicates",
#          "180 Peptides | 9 Peptides Per Pool | 3 Replicates",
#          "360 Peptides | 9 Peptides Per Pool | 3 Replicates",
#          "720 Peptides | 9 Peptides Per Pool | 3 Replicates")
# )

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
    # geom_hline(yintercept = num.first.round.pools, linetype = "dashed") +
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
         filename = paste0(OUTPUT.DIR, '/', num.peptides, 'peptides_', num.peptides.per.pool, 'perpool_', num.coverage, 'x.pdf'),
         width = 9, height = 6, dpi = 300
  )
}
