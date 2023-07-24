library(ggplot2)
library(ggsignif)


MERGED.TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace_vs_naive_approach/deconvolutethis_designs/ace_vs_naive_approach_benchmark_experiment_deconvolutethis_designs_results_merged.tsv"
OUTPUT.DIR <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace_vs_naive_approach/deconvolutethis_designs"

df.plot <- read.csv(MERGED.TSV.FILE, sep = '\t')
df.plot$group <- paste0(df.plot$num_peptides, "/", df.plot$num_peptides_per_pool, "/", df.plot$num_coverage)
df.plot$num_positive_peptide_sequences <- as.character(df.plot$num_positive_peptide_sequences)

# Stratified
plot <- ggplot(df.plot, aes(x = num_positive_peptide_sequences, y = predicted_total_pools, fill = solver)) +
  geom_boxplot(width = 0.382, outlier.size = 0.382, notchwidth = 0.382) +
  facet_wrap(~ group, scales = "free_y", ncol = 1) +
  scale_x_discrete(limits = c("1","2","3","4","5","6","7","8","9","10","15","20","25","30","40"))
ggsave(plot = plot, filename = paste0(OUTPUT.DIR, "/ace_vs_naive_approach_deconvolute_designs_results.pdf"),
       width = 32, height = 96, dpi = 300, limitsize = F,
)
# geom_signif(comparisons = list(c("ace_1", "ace_2", "bogey")), 
#             map_signif_level = TRUE) +

# 
# plot <- ggplot(df.plot, aes(x = solver, y = num_total_pools,)) +
#   geom_boxplot(width = 0.382, outlier.size = 0.382, notchwidth = 0.382) +
#   scale_y_continuous(limits = c(0,180)) +
#   facet_wrap(~ group, scales = "free_y", ncol = 2)
# print(plot)
