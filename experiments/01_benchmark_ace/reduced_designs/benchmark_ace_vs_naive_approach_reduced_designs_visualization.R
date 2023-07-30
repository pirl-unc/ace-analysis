library(ggplot2)
library(ggsignif)


DATA.DIR <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace_vs_naive_approach/reduced_designs"

df.plot <- data.frame()
for (tsv.file in Sys.glob(paste0(DATA.DIR, "/*/*results*.tsv"))) {
  df.temp <- read.csv(tsv.file, sep = "\t")
  df.plot <- rbind(df.plot, df.temp)
}

df.plot$group <- paste0(df.plot$num_peptides, "/", df.plot$num_peptides_per_pool, "/", df.plot$num_coverage)
df.plot$num_positive_peptide_sequences <- as.character(df.plot$num_positive_peptide_sequences)

# Stratified
plot <- ggplot(df.plot, aes(x = num_positive_peptide_sequences, y = predicted_total_pools, fill = solver)) +
  geom_boxplot(width = 0.382, outlier.size = 0.382, notchwidth = 0.382) +
  facet_wrap(~ group, scales = "free", ncol = 1) +
  scale_x_discrete(limits = c(
    "1", "2", "4", "8", "16", "32"
  ))
print(plot)
