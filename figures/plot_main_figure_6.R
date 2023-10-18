library(ggplot2)
library(ggpubr)
library(dplyr)


# Step 1. Define constants
OUTPUT.DIR <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/figures"
POINTRANGE.COLOR.HEX <- "#000000"
BAR.COLOR.HEX <- "#DCDCDC"
BAR.WIDTH <- 0.5
POINTRANGE.SIZE <- 1.0
COMMON.THEME <- theme(axis.title = element_text(size = 14),
                      axis.text = element_text(size = 12),
                      legend.title = element_text(size = 12),
                      legend.text = element_text(size = 12))

# Step 2. Plot panel A
TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/02_benchmark_ace_runtime_and_memory/runtime_memory_benchmark_experiment_results_merged.tsv"
df <- read.csv(TSV.FILE, sep = "\t")
df.plot <- df %>% 
  dplyr::filter(num_peptides_per_pool == 10) %>%
  dplyr::filter(num_coverage == 3)
df.plot$num_peptides <- as.character(df.plot$num_peptides)
df.plot$num_peptides <- factor(
  df.plot$num_peptides,
  levels = c("100","200","300","400","500",
             "600","700","800","900","1000")
)
figure.panel.a <- ggplot(df.plot, aes(x = num_peptides, y = runtime_seconds)) +
  geom_pointrange(stat = "summary", 
                  fun.min = min, 
                  fun.max = max, 
                  fun = median, 
                  size = POINTRANGE.SIZE, 
                  color = POINTRANGE.COLOR.HEX) + 
  xlab("Number of Peptides") + ylab("Runtime (Seconds)") + ggtitle("10 Peptides per Pool | 3x Coverage") +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0,50,5)) +
  theme_pubr() + COMMON.THEME

# Step 3. Plot panel B
df.plot.temp <- df.plot %>%
  dplyr::group_by(num_peptides) %>%
  dplyr::summarise(peak_memory_mb = mean(peak_memory_mb))
figure.panel.b <- ggplot(df.plot.temp, aes(x = num_peptides, y = peak_memory_mb)) +
  geom_col(fill = BAR.COLOR.HEX, color = "black", width = BAR.WIDTH) +
  xlab("Number of Peptides") + ylab("Average Peak Memory (Megabytes)") + ggtitle("10 Peptides per Pool | 3x Coverage") +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0,8,2)) +
  theme_pubr() + COMMON.THEME

# Step 4. Plot panel C
df.plot <- df %>% 
  dplyr::filter(num_peptides == 1000) %>%
  dplyr::filter(num_coverage == 3)
df.plot$num_peptides_per_pool <- as.character(df.plot$num_peptides_per_pool)
df.plot$num_peptides_per_pool <- factor(
  df.plot$num_peptides_per_pool,
  levels = c("10","12","14","16","18",
             "20","22","24","26","28","30")
)
figure.panel.c <- ggplot(df.plot, aes(x = num_peptides_per_pool, y = runtime_seconds)) +
  geom_pointrange(stat = "summary", 
                  fun.min = min, 
                  fun.max = max, 
                  fun = median, 
                  size = POINTRANGE.SIZE,
                  color = POINTRANGE.COLOR.HEX) + 
  xlab("Number of Peptides per Pool") + ylab("Runtime (Seconds)") + ggtitle("1000 Peptides | 3x Coverage") +
  scale_y_continuous(limits = c(0, 70), breaks = seq(0,70,10)) +
  theme_pubr() + COMMON.THEME

# Step 5. Plot panel D
df.plot.temp <- df.plot %>%
  dplyr::group_by(num_peptides_per_pool) %>%
  dplyr::summarise(peak_memory_mb = mean(peak_memory_mb))
figure.panel.d <- ggplot(df.plot.temp, aes(x = num_peptides_per_pool, y = peak_memory_mb)) +
  geom_col(fill = BAR.COLOR.HEX, color = "black", width = BAR.WIDTH) +
  xlab("Number of Peptides per Pool") + ylab("Average Peak Memory (Megabytes)") + ggtitle("1000 Peptides | 3x Coverage") +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0,20,5)) +
  theme_pubr() + COMMON.THEME

# Step 6. Plot panel E
df.plot <- df %>% 
  dplyr::filter(num_peptides == 1000) %>%
  dplyr::filter(num_peptides_per_pool == 10)
df.plot$num_coverage <- as.character(df.plot$num_coverage)
df.plot$num_coverage <- factor(
  df.plot$num_coverage,
  levels = c("3","4","5","6","7","8","9","10")
)
figure.panel.e <- ggplot(df.plot, aes(x = num_coverage, y = runtime_seconds)) +
  geom_pointrange(stat = "summary", 
                  fun.min = min, 
                  fun.max = max, 
                  fun = median, 
                  size = POINTRANGE.SIZE,
                  color = POINTRANGE.COLOR.HEX) + 
  xlab("Coverage") + ylab("Runtime (Seconds)") + ggtitle("1000 Peptides | 10 Peptides per Pool") +
  scale_y_continuous(limits = c(0, 70), breaks = seq(0,70,5)) +
  theme_pubr() + COMMON.THEME

# Step 7. Plot panel F
df.plot.temp <- df.plot %>%
  dplyr::group_by(num_coverage) %>%
  dplyr::summarise(peak_memory_mb = mean(peak_memory_mb))
figure.panel.f <- ggplot(df.plot.temp, aes(x = num_coverage, y = peak_memory_mb)) +
  geom_col(fill = BAR.COLOR.HEX, color = "black", width = BAR.WIDTH) +
  xlab("Coverage") + ylab("Average Peak Memory (Megabytes)") + ggtitle("1000 Peptides | 10 Peptides per Pool") +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20)) +
  theme_pubr() + COMMON.THEME

# Step 8. Plot figure 6
figure.6 <- ggarrange(plotlist = list(figure.panel.a,
                                      figure.panel.c,
                                      figure.panel.e,
                                      figure.panel.b,
                                      figure.panel.d,
                                      figure.panel.f),
                      ncol = 3, nrow = 2, align = "hv",
                      widths = c(1,1))
ggsave(plot = figure.6, filename = paste0(OUTPUT.DIR, "/main_figure_6.pdf"),
       width = 16, height = 9, dpi = 300)
