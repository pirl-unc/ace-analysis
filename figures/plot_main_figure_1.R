library(ggplot2)
library(dplyr)
library(ggpubr)


# Step 1. Define constants
OUTPUT.DIR <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/figures"
COMMON.THEME <- theme(axis.title = element_text(size = 14),
                      axis.text = element_text(size = 12),
                      legend.title = element_text(size = 12),
                      legend.text = element_text(size = 12))

# Step 2. Plot panel (a)
CSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/precision_recall_violations_iterations/precision_recall_f1_violations.csv"
df.plot <- read.csv(CSV.FILE)
plot.precision <- ggplot(df.plot, aes(x = iteration, y = precision)) +
  geom_line() +
  xlab("Iteration") + ylab("Precision") +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1.0)) +
  theme_pubr() + COMMON.THEME
plot.recall <- ggplot(df.plot, aes(x = iteration, y = recall)) +
  geom_line() +
  xlab("Iteration") + ylab("Sensitivity") +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1.0)) +
  theme_pubr() + COMMON.THEME
plot.num.violations <- ggplot(df.plot, aes(x = iteration, y = num_violations)) +
  geom_line() +
  xlab("Iteration") + ylab("Number of Violations") +
  scale_y_continuous(limits = c(0,14000), breaks = c(0,3500,7000,10500,14000)) +
  theme_pubr() + COMMON.THEME
figure <- ggarrange(plotlist = list(plot.precision, plot.recall, plot.num.violations),
                    ncol = 1, align = "hv")
print(figure)

# Step 3. Plot panel (b)
MERGED.TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/reduced_designs/reduced_designs_experiment_results_merged.tsv"


ggsave(filename = paste0(OUTPUT.DIR, "/precision_recall_violations_iterations.pdf"),
       width = 8, height = 8, dpi = 300)


