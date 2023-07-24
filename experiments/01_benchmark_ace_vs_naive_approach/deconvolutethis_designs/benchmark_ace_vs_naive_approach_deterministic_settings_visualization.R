library(ggplot2)
library(ggsignif)


BENCHMARK.EXPERIMENTS.CSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace_vs_naive_approach/ace_vs_naive_approach_benchmark_experiment_results.csv"

df.plot <- read.csv(BENCHMARK.EXPERIMENTS.CSV.FILE)
df.plot$group <- paste0(df.plot$num_peptides, "/", df.plot$num_peptides_per_pool, "/", df.plot$num_coverage)
df.plot$num_true_positive_peptides <- as.character(df.plot$num_true_positive_peptides)

# Stratified
plot <- ggplot(df.plot, aes(x = num_true_positive_peptides, y = num_total_pools, fill = solver)) +
  geom_boxplot(width = 0.382, outlier.size = 0.382, notchwidth = 0.382) +
  facet_wrap(~ group, scales = "free_y", ncol = 2) +
  scale_x_discrete(limits = c("1","2","3","4","5","6","7","8","9","10","15","20","25","30","40"))
print(plot)
# geom_signif(comparisons = list(c("ace_1", "ace_2", "bogey")), 
#             map_signif_level = TRUE) +


plot <- ggplot(df.plot, aes(x = solver, y = num_total_pools,)) +
  geom_boxplot(width = 0.382, outlier.size = 0.382, notchwidth = 0.382) +
  scale_y_continuous(limits = c(0,180)) +
  facet_wrap(~ group, scales = "free_y", ncol = 2)
print(plot)
