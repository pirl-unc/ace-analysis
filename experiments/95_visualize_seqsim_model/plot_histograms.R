library(ggplot2)
library(ggpubr)


CSV.FILE.1 <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/seqsim_model/outgroup_euclidean_not_finetuned.csv"
CSV.FILE.2 <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/seqsim_model/ingroup_euclidean_not_finetuned.csv"
CSV.FILE.3 <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/seqsim_model/outgroup_euclidean_finetuned.csv"
CSV.FILE.4 <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/seqsim_model/ingroup_euclidean_finetuned.csv"
OUTPUT.DIR <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/seqsim_model"
OUTGROUP.HEX <- "#173ccc"
INGROUP.HEX <- "#e1151b"
FILL.COLORS <- c("ingroup" = INGROUP.HEX,
                 "outgroup" = OUTGROUP.HEX)

df.outgroup.not.finetuned <- read.csv(CSV.FILE.1)
df.outgroup.not.finetuned$group <- "outgroup"
df.ingroup.not.finetuned <- read.csv(CSV.FILE.2)
df.ingroup.not.finetuned$group <- "ingroup"
df.not.finetuned <- rbind(df.outgroup.not.finetuned, df.ingroup.not.finetuned)

df.outgroup.finetuned <- read.csv(CSV.FILE.3)
df.outgroup.finetuned$group <- "outgroup"
df.ingroup.finetuned <- read.csv(CSV.FILE.4)
df.ingroup.finetuned$group <- "ingroup"
df.finetuned <- rbind(df.outgroup.finetuned, df.ingroup.finetuned)

plot.not.finedtuned <- ggplot(df.not.finetuned, aes(x = X0, fill = group)) +
  geom_density(alpha = 0.618) +
  xlab("Euclidean Similarity") + ylab("Density") +
  scale_x_continuous(limits = c(0,1)) +
  scale_fill_manual(values = FILL.COLORS) +
  theme_pubr()
print(plot.not.finedtuned)
ggsave(plot = plot.not.finedtuned, filename = paste0(OUTPUT.DIR, "/density_plot_not_fine-tuned.pdf"),
       width = 16, height = 8, dpi = 300)

plot.finedtuned <- ggplot(df.finetuned, aes(x = X0, fill = group)) +
  geom_density(alpha = 0.618) +
  xlab("Euclidean Similarity") + ylab("Density") +
  scale_x_continuous(limits = c(0,1)) +
  scale_fill_manual(values = FILL.COLORS) +
  theme_pubr()
print(plot.finedtuned)
ggsave(plot = plot.finedtuned, filename = paste0(OUTPUT.DIR, "/density_plot_fine-tuned.pdf"),
       width = 16, height = 8, dpi = 300)
