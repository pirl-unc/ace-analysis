library(ggplot2)


TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/scripts/ace-analysis/examples/test_deconvolution_10percent_randomeffectsoff_results.tsv"

df <- read.csv(TSV.FILE, sep = '\t')
df.plot <- data.frame(
  auroc = c(df$aucroc_score_empirical, df$aucroc_score_lasso, df$aucroc_score_em),
  method = c(rep('empirical', nrow(df)), rep('lasso', nrow(df)), rep('em', nrow(df)))
)
plot <- ggplot(df.plot, aes(x = auroc, fill = method)) +
  geom_histogram(position="identity", binwidth = 0.01, alpha = 0.5) +
  geom_density(alpha = 0.6)
print(plot)


plot <- ggplot(df.plot, aes(x = method, y = auroc)) +
  geom_boxplot(position="identity")
print(plot)
