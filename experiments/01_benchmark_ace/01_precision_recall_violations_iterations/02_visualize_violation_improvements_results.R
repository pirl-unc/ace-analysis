library(ggplot2)
library(ggpubr)
library(dplyr)


# Step 1. Define constants
CSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/01_precision_recall_violations_iterations/precision_recall_violations.csv"
OUTPUT.DIR <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/01_precision_recall_violations_iterations"
df <- read.csv(CSV.FILE)
df.plot <- df %>%
  group_by(iteration) %>%
  summarise(precision.mean = mean(precision),
            precision.std = sd(precision),
            recall.mean = mean(recall),
            recall.std = sd(recall),
            f1.mean = mean(f1),
            f1.std = sd(f1),
            violations.mean = mean(num_violations),
            violations.std = sd(num_violations))

df.plot.precision <- df.plot %>%
  dplyr::reframe(iteration = iteration,
                 precision.mean = precision.mean,
                 precision.lower = ifelse(precision.mean - precision.std < 0, 0, precision.mean - precision.std),
                 precision.upper = ifelse(precision.mean + precision.std > 1, 1, precision.mean + precision.std))
plot.precision <- ggplot(df.plot.precision, aes(x = iteration)) +
  geom_ribbon(aes(ymax = precision.upper, ymin = precision.lower), fill = "#fcf3f8") +
  geom_line(aes(y = precision.mean), color = "#a22867", linewidth = 0.618) +
  geom_line(aes(y = precision.lower), color = "#d1458e", linewidth = 0.618) +
  geom_line(aes(y = precision.upper), color = "#e08284", linewidth = 0.618) +
  xlab("Iteration") + ylab("Precision") +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1.0)) +
  theme_pubr()
print(plot.precision)

df.plot.recall <- df.plot %>%
  dplyr::reframe(iteration = iteration,
                 recall.mean = recall.mean,
                 recall.lower = ifelse(recall.mean - recall.std < 0, 0, recall.mean - recall.std),
                 recall.upper = ifelse(recall.mean + recall.std > 1, 1, recall.mean + recall.std))
plot.recall <- ggplot(df.plot.recall, aes(x = iteration)) +
  geom_ribbon(aes(ymax = recall.upper, ymin = recall.lower), fill = "#fcf3f8") +
  geom_line(aes(y = recall.mean), color = "#a22867", linewidth = 0.618) +
  geom_line(aes(y = recall.lower), color = "#d1458e", linewidth = 0.618) +
  geom_line(aes(y = recall.upper), color = "#e08284", linewidth = 0.618) +
  xlab("Iteration") + ylab("Sensitivity") +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1.0)) +
  theme_pubr()
print(plot.recall)

df.plot.violations <- df.plot %>%
  dplyr::reframe(iteration = iteration,
                 violations.mean = violations.mean,
                 violations.lower = ifelse(violations.mean - violations.std < 0, 0, violations.mean - violations.std),
                 violations.upper = ifelse(violations.mean + violations.std > 1, 1, violations.mean + violations.std))
plot.violations <- ggplot(df.plot.violations, aes(x = iteration)) +
  geom_line(aes(y = violations.mean), color = "#a22867", linewidth = 0.618) +
  xlab("Iteration") + ylab("Number of Violations") +
  scale_y_continuous(limits = c(0,8000), breaks = c(0,2000,4000,6000,8000)) +
  theme_pubr()
print(plot.violations)

figure <- ggarrange(plotlist = list(plot.precision, plot.recall, plot.violations),
                    ncol = 1, align = "hv")
print(figure)
ggsave(filename = paste0(OUTPUT.DIR, "/precision_recall_violations_iterations.pdf"),
       width = 8, height = 8, dpi = 300)
