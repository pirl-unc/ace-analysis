library(ggplot2)
library(ggpubr)
library(dplyr)
library(viridis)


# Step 1. Define constants
OUTPUT.DIR <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/figures"
COMMON.THEME <- theme(axis.title = element_text(size = 14),
                      axis.text = element_text(size = 12),
                      legend.title = element_text(size = 12),
                      legend.text = element_text(size = 12))
PRECISION.LOW.HEX <- "#eaecee"
PRECISION.HIGH.HEX <- "#2f4858"
NUM.TOTAL.POOLS.LOW.HEX <- "#e5f2ef"
NUM.TOTAL.POOLS.HIGH.HEX <- "#008563"
LWD <- 0.618
LINETYPE <- 1
VALUE.FONT.SIZE <- 3.6
NUM.PEPTIDES.PER.POOL <- c(5,8,10,20,25,40)
TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/07_design_exploration/design_exploration_experiment_results_merged.tsv"
df <- read.csv(TSV.FILE, sep = "\t")

PlotPrecision <- function(df.plot, title) {
  df.plot$num_peptides_per_pool <- as.character(df.plot$num_peptides_per_pool)
  df.plot$num_peptides_per_pool <- factor(
    df.plot$num_peptides_per_pool,
    levels = as.character(NUM.PEPTIDES.PER.POOL)
  )
  df.plot$text_color <- ifelse(df.plot$precision < 0.5, "white", "black")
  plot <- ggplot(df.plot, aes(x = num_peptides_per_pool, y = num_coverage, fill = precision)) +
    geom_tile(color = "white", lwd = LWD, linetype = LINETYPE) +
    geom_text(aes(label = precision, color = text_color), size = VALUE.FONT.SIZE) +
    xlab("Number of Peptides per Pool") + ylab("Coverage") + ggtitle(title) +
    scale_fill_viridis(option = "magma", direction = 1, limits = c(0,1)) +
    scale_x_discrete() +
    scale_y_continuous(breaks = seq(2,10,1)) +
    scale_color_manual(values = c('black', 'white'))+
    theme_pubr() + COMMON.THEME
  return(plot)
}

PlotNumTotalPools <- function(df.plot, title, min, max) {
  df.plot$text_color <- ifelse(df.plot$num_total_pools < mean(c(min(df.plot$num_total_pools), max(df.plot$num_total_pools))), "black", "white")
  df.plot$num_peptides_per_pool <- factor(
    df.plot$num_peptides_per_pool,
    levels = as.character(NUM.PEPTIDES.PER.POOL)
  )
  plot <- ggplot(df.plot, aes(x = num_peptides_per_pool, y = num_coverage, fill = num_total_pools)) +
    geom_tile(color = "white", lwd = LWD, linetype = LINETYPE) +
    geom_text(aes(label = num_total_pools, color = text_color), size = VALUE.FONT.SIZE) +
    xlab("Number of Peptides per Pool") + ylab("Coverage") + ggtitle(title) +
    scale_fill_viridis(option = "viridis", direction = -1, limits = c(min, max)) +
    scale_x_discrete() +
    scale_y_continuous(breaks = seq(2,10,1)) +
    scale_color_manual(values = c('black', 'white'))+
    theme_pubr() + COMMON.THEME
  return(plot)
}


# Step 2. Plot panel A
df.plot <- df %>%
  filter(num_peptides == 200) %>%
  filter(num_peptides_per_pool %in% NUM.PEPTIDES.PER.POOL) %>%
  group_by(num_peptides_per_pool, num_coverage) %>%
  summarise(precision = mean(precision_empirical),
            num_total_pools = mean(predicted_total_pools))
df.plot$precision <- round(df.plot$precision, 2)
df.plot$num_total_pools <- round(df.plot$num_total_pools)
panel.a.precision <- PlotPrecision(df.plot = df.plot, title = "200 Peptides")
panel.a.pools <- PlotNumTotalPools(df.plot = df.plot, title = "200 Peptides", min = 100, max = 500)


# Step 3. Plot panel B
df.plot <- df %>%
  filter(num_peptides == 400) %>%
  filter(num_peptides_per_pool %in% NUM.PEPTIDES.PER.POOL) %>%
  group_by(num_peptides_per_pool, num_coverage) %>%
  summarise(precision = mean(precision_empirical),
            num_total_pools = mean(predicted_total_pools))
df.plot$precision <- round(df.plot$precision, 2)
df.plot$num_total_pools <- round(df.plot$num_total_pools)
panel.b.precision <- PlotPrecision(df.plot = df.plot, title = "400 Peptides")
panel.b.pools <- PlotNumTotalPools(df.plot = df.plot, title = "400 Peptides", min = 200, max = 900)


# Step 4. Plot panel C
df.plot <- df %>%
  filter(num_peptides == 600) %>%
  filter(num_peptides_per_pool %in% NUM.PEPTIDES.PER.POOL) %>%
  group_by(num_peptides_per_pool, num_coverage) %>%
  summarise(precision = mean(precision_empirical),
            num_total_pools = mean(predicted_total_pools))
df.plot$precision <- round(df.plot$precision, 2)
df.plot$num_total_pools <- round(df.plot$num_total_pools)
panel.c.precision <- PlotPrecision(df.plot = df.plot, title = "600 Peptides")
panel.c.pools <- PlotNumTotalPools(df.plot = df.plot, title = "600 Peptides", min = 300, max = 1300)


# Step 5. Plot panel D
df.plot <- df %>%
  filter(num_peptides == 800) %>%
  filter(num_peptides_per_pool %in% NUM.PEPTIDES.PER.POOL) %>%
  group_by(num_peptides_per_pool, num_coverage) %>%
  summarise(precision = mean(precision_empirical),
            num_total_pools = mean(predicted_total_pools))
df.plot$precision <- round(df.plot$precision, 2)
df.plot$num_total_pools <- round(df.plot$num_total_pools)
panel.d.precision <- PlotPrecision(df.plot = df.plot, title = "800 Peptides")
panel.d.pools <- PlotNumTotalPools(df.plot = df.plot, title = "800 Peptides", min = 400, max = 1700)


# Step 6. Plot figure
figure <- ggarrange(plotlist = list(panel.a.precision,
                                    panel.a.pools,
                                    panel.b.precision,
                                    panel.b.pools,
                                    panel.c.precision,
                                    panel.c.pools,
                                    panel.d.precision,
                                    panel.d.pools),
                    ncol = 2, nrow = 4, align = "hv", common.legend = F)
ggsave(filename = paste0(OUTPUT.DIR, "/main_figure_5.pdf"), width = 8, height = 16, dpi = 300)
