library(ggplot2)
library(ggsignif)


TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace_vs_naive_approach/alanine_scanning/ace_vs_naive_approach_alanine_scanning_benchmark_experiment_results_merged.tsv"

df.plot <- read.csv(TSV.FILE, sep = '\t')
df.plot$group <- paste0(df.plot$num_peptides, "/", df.plot$num_peptides_per_pool, "/", df.plot$num_coverage)
df.plot$num_positive_peptide_sequences <- as.character(df.plot$num_positive_peptide_sequences)

# Stratified
plot <- ggplot(df.plot, aes(x = num_positive_peptide_sequences, y = predicted_total_pools, fill = solver)) +
  geom_boxplot(width = 0.382, outlier.size = 0.382, notchwidth = 0.382) +
  facet_wrap(~ group, scales = "free_y", ncol = 3) +
  scale_x_discrete(limits = c(
    "9", "18", "27", "36", "45", "54",
    "63", "72", "81", "90", "99", "108",
    "117", "126", "135", "144", "153", "162",
    "171", "180", "189", "198", "207", "216"
  ))
print(plot)
# geom_signif(comparisons = list(c("ace_1", "ace_2", "bogey")), 
#             map_signif_level = TRUE) +


plot <- ggplot(df.plot, aes(x = solver, y = num_total_pools,)) +
  geom_boxplot(width = 0.382, outlier.size = 0.382, notchwidth = 0.382) +
  scale_y_continuous(limits = c(0,180)) +
  facet_wrap(~ group, scales = "free_y", ncol = 2)
print(plot)
