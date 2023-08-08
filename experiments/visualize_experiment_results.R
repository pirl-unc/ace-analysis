library(ggplot2)
library(plyr)


MERGED.TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/scripts/ace-analysis/examples/177peptides_10perpool_3x/tan_et_al_inf_gen_evol_2021_experiment_results_295153904.tsv"
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

plot.aucroc.score <- ggplot(df.plot, aes(x = num_positive_peptide_sequences,
                                        y = precision_empirical,
                                        fill = solver)) +
  xlab("Percentage of Positive Peptides") + ylab("AUCROC Score") +
  geom_boxplot() +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12))
print(plot.aucroc.score)

plot.specificity <- ggplot(df.plot, aes(x = num_positive_peptide_sequences,
                                        y = specificity_empirical,
                                        fill = solver)) +
  xlab("Percentage of Positive Peptides") + ylab("AUCROC Score") +
  geom_boxplot() +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12))
print(plot.specificity)

plot.precision <- ggplot(df.plot, aes(x = num_positive_peptide_sequences,
                                      y = precision_empirical,
                                      fill = solver)) +
  xlab("Percentage of Positive Peptides") + ylab("Precision") +
  geom_boxplot() +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12))
print(plot.precision)


plot.num.violations <- ggplot(df.plot, aes(x = num_positive_peptide_sequences,
                                           y = num_violations,
                                           fill = solver)) +
  xlab("AUROC Score") + ylab("Num Violations") +
  geom_boxplot() +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12))
print(plot.num.violations)
