library(ggplot2)
library(ggpubr)
library(dplyr)


CSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/precision_recall_violations_iterations/precision_recall_f1_violations.csv"
OUTPUT.DIR <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/precision_recall_violations_iterations"
df.plot <- read.csv(CSV.FILE)

# Find mean for each iteration
df.plot <- df.plot %>%
  group_by(iteration) %>%
  summarise(precision.mean = mean(precision))
p <- ggplot(df.plot, aes(x = iteration, y = precision.mean)) +
  geom_line()
print(p)

plot.precision <- ggplot(df.plot, aes(x = iteration, y = precision)) +
  geom_line() +
  xlab("Iteration") + ylab("Precision") +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1.0)) +
  theme_pubr()
plot.recall <- ggplot(df.plot, aes(x = iteration, y = recall)) +
  geom_line() +
  xlab("Iteration") + ylab("Sensitivity") +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1.0)) +
  theme_pubr()
plot.num.violations <- ggplot(df.plot, aes(x = iteration, y = num_violations)) +
  geom_line() +
  xlab("Iteration") + ylab("Number of Violations") +
  scale_y_continuous(limits = c(0,14000), breaks = c(0,3500,7000,10500,14000)) +
  theme_pubr()
figure <- ggarrange(plotlist = list(plot.precision, plot.recall, plot.num.violations),
                    ncol = 1, align = "hv")
print(figure)
ggsave(filename = paste0(OUTPUT.DIR, "/precision_recall_violations_iterations.pdf"),
       width = 8, height = 8, dpi = 300)
