library(ggplot2)
library(plyr)


MERGED.TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/scripts/ace-analysis/examples/test_results.tsv"
df.plot <- read.csv(MERGED.TSV.FILE, sep = '\t')
df.plot$num_positive_peptide_sequences <- as.character(df.plot$num_positive_peptide_sequences)
plot.num.total.pools <- ggplot(df.plot, aes(x = num_positive_peptide_sequences,
                                            y = predicted_total_pools,
                                            fill = solver)) +
  xlab("Percentage of Positive Peptides") + ylab("Number of Total Pools") +
  geom_boxplot() +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12))
print(plot.num.total.pools)

