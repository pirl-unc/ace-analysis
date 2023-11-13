library(ggplot2)
library(plyr)
library(dplyr)
library(ggpubr)


# Step 1. Define constants
OUTPUT.DIR <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/figures"
ACE.HEX <- "#d65db1"
DECONVOLUTETHIS.HEX <- "#53424C"
COMMON.THEME <- theme(axis.title = element_text(size = 14),
                      axis.text = element_text(size = 12),
                      legend.title = element_text(size = 14),
                      legend.text = element_text(size = 12))


# Step 2. Plot
MERGED.TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/04_deconvolutethis_designs/deconvolutethis_experiment_results_merged.tsv"
DECONVOLUTETHIS.DESIGNS.CSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/raw/deconvolutethis_designs.csv"
df <- read.csv(MERGED.TSV.FILE, sep = '\t')
df.deconvolutethis <- read.csv(DECONVOLUTETHIS.DESIGNS.CSV.FILE)
df.deconvolutethis$group <- paste0(df.deconvolutethis$num_peptides, "/", df.deconvolutethis$num_peptides_per_pool, "/", df.deconvolutethis$num_coverage)
df$group <- paste0(df$num_peptides, "/", df$num_peptides_per_pool, "/", df$num_coverage)
df.deconvolutethis$num_positive_peptide_sequences <- df.deconvolutethis$num_true_positive_peptides
df.plot <- df %>%
  dplyr::group_by(group, num_positive_peptide_sequences) %>%
  dplyr::summarise(num_total_pools = round(mean(predicted_total_pools)))
df.plot$method <- "ACE"
df.plot <- rbind(df.plot, data.frame(
  group = df.deconvolutethis$group,
  num_positive_peptide_sequences = df.deconvolutethis$num_positive_peptide_sequences,
  num_total_pools = df.deconvolutethis$num_deconvolutethis_total_pools,
  method = "DeconvoluteThis"
))
df.plot$num_positive_peptide_sequences <- as.character(df.plot$num_positive_peptide_sequences)
df.plot$num_positive_peptide_sequences <- factor(
  df.plot$num_positive_peptide_sequences,
  levels = c(
    "1","2","3","4","5","6","7","8","9","10","15","20","25","30","40"
  )
)
df.plot$group <- factor(
  df.plot$group,
  levels = c(
    "120/12/3",
    "120/15/3",
    "120/20/3",
    "120/24/3",
    "800/16/4",
    "800/20/3",
    "800/20/4",
    "800/32/4",
    "800/32/5",
    "800/50/4",
    "800/50/5",
    "800/80/5",
    "800/80/6",
    "800/100/4",
    "800/100/5",
    "800/100/6",
    "800/160/4",
    "800/160/5"
  )
)
barplot <- ggplot(df.plot, aes(x = num_positive_peptide_sequences, y = num_total_pools, fill = method)) +
  geom_bar(position="dodge", stat="identity", color = "black", width = 0.618) +
  xlab("Number of Positive Peptides") + ylab("Total Number of Pools") +
  guides(fill = guide_legend(title = "Method")) +
  facet_wrap(~group, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = c("ACE" = ACE.HEX, "DeconvoluteThis" = DECONVOLUTETHIS.HEX)) +
  theme_bw() + COMMON.THEME +
  theme(panel.grid.minor.y = element_blank(),
        strip.text = element_text(size = 12))
print(barplot)
ggsave(plot = barplot, filename = paste0(OUTPUT.DIR, "/supplementary_figure_1.pdf"),
       width = 16, height = 16, dpi = 300)
