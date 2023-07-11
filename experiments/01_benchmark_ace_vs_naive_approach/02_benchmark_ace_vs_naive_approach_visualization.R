library(ggplot2)
library(ggsignif)


BENCHMARK.EXPERIMENTS.CSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/03_benchmark_ace_vs_naive_approach/ace_vs_naive_approach_benchmark_experiment_results.csv"

df.plot <- read.csv(BENCHMARK.EXPERIMENTS.CSV.FILE)

plot <- ggplot(df.plot, aes(x = solver, y = num_total_pools, fill = solver)) +
  geom_boxplot() +
  geom_signif(comparisons = list(c("ace", "bogey")), 
              map_signif_level=TRUE) +
  geom_hline(yintercept = 67, linetype = "dashed")

print(plot)
