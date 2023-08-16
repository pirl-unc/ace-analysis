library(ggplot2)
library(plyr)
library(dplyr)
library(ggpubr)
library(ggsignif)


MERGED.TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/real_world_datasets/cameron_et_al_sci_trans_med_2013//cameron2013_dataset_experiment_results_merged.tsv"
OUTPUT.DIR <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/real_world_datasets/cameron_et_al_sci_trans_med_2013/"
DODGE.WIDTH <- 0.5
ERRORBAR.WIDTH <- 0.146
BOXPLOT.WIDTH <- 0.263
FILL.COLORS <- c("ACE" = "#845EC2",
                 "ACE (w/o Clustering)" = "#D65DB1",
                 "Randomized Block Design" = "#FF9671", # F05B61
                 "Repeated Block Design" = "#ABABAA")

# Step 2. Load data
df.plot <- read.csv(MERGED.TSV.FILE, sep = '\t')

# Step 3. Postprocess data
df.plot$group <- paste0(df.plot$num_peptides, "/", df.plot$num_peptides_per_pool, "/", df.plot$num_coverage)
df.plot <- df.plot[df.plot$solver %in% c("ace_golfy_clusteron_noextrapools",
                                         "ace_golfy_clusteroff_noextrapools",
                                         "randomized_block_assignment",
                                         "repeated_block_assignment"),]

# Step 4. Map values
df.plot$solver <- mapvalues(
  x = df.plot$solver, 
  from = c("ace_golfy_clusteron_noextrapools",
           "ace_golfy_clusteroff_noextrapools",
           "randomized_block_assignment",
           "repeated_block_assignment"),
  to = c("ACE",
         "ACE (w/o Clustering)",
         "Randomized Block Design",
         "Repeated Block Design")
)
df.plot$solver <- factor(
  df.plot$solver,
  levels = c("ACE",
             "ACE (w/o Clustering)",
             "Randomized Block Design",
             "Repeated Block Design")
)

# Step 5. Plot
num.peptides <- df.plot$num_peptides[1]
num.peptides.per.pool <- df.plot$num_peptides_per_pool[1]
num.coverage <- df.plot$num_coverage[1]
num.first.round.pools <- ceiling(df.plot$num_peptides[1] / df.plot$num_peptides_per_pool[1]) * df.plot$num_coverage[1]
num.max.pools <- max(df.plot$predicted_total_pools)
plot.num.total.pools <- ggplot(df.plot, aes(x = solver, y = predicted_total_pools, fill = solver)) + 
  geom_hline(yintercept = num.first.round.pools, linetype = "dashed") +
  stat_boxplot(geom = "errorbar", width = 0.4) +
  geom_boxplot(outlier.shape = NA, width = 0.8) +
  xlab("") + ylab("Number of Total Pools") + 
  scale_fill_manual(values = FILL.COLORS) +
  scale_y_continuous(limits = c(18, 36), breaks = c(18,27,36)) +
  theme_pubr()
ggsave(plot = plot.num.total.pools, filename = paste0(OUTPUT.DIR, "/cameron2013_total_pools.pdf"),
       width = 4, height = 8, dpi = 300)

df.plot.temp <- df.plot %>%
  dplyr::group_by(solver) %>%
  dplyr::summarise(
    sd = sd(precision_empirical),
    precision = mean(precision_empirical)
  )
plot.precision <- ggplot(df.plot.temp, aes(x = solver, y = precision, fill = solver)) +
  geom_col(width = 0.8, colour = "black", position = "dodge") +
  geom_errorbar(aes(ymin = precision - sd, ymax = ifelse(precision + sd > 1, 1, precision + sd)), width = 0.2, position = position_dodge(0.5)) +
  xlab("") + ylab("Precision") +
  scale_y_continuous(limits = c(0,1), breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) +
  scale_fill_manual(values = FILL.COLORS) +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  theme_pubr()
print(plot.precision)
ggsave(plot = plot.precision, filename = paste0(OUTPUT.DIR, "/cameron2013_precision.pdf"),
       width = 4, height = 8, dpi = 300)



