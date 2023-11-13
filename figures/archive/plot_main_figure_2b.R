library(ggplot2)
library(dplyr)
library(plyr)
library(ggpubr)
library(ggsignif)


EVALUATION.RESULTS.DIR <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/02_benchmark_ace/01_benchmark_against_existing_tools/06_evaluation_results/300immunogenic_5nonimmunogenic_0dispersion_fnr0.00"
DECONVOLUTETHIS.DESIGN.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/raw/references/deconvolutethis_designs.csv"
ERRORBAR.WIDTH <- 0.384
BOXPLOT.WIDTH <- 0.618
# SOLVER.COLORS <- c("ace-ace-em" = "#A500B1",
#                    "randomized-ace-em" = "#FF9671",
#                    "repeated-ace-em" = "#ABABAA",
#                    "deconvolutethis-deconvolutethis" = "#AA85AA",
#                    "rickard-rickard-background-subtracted" = "#6b76eb",
#                    "rickard-rickard_bayesian" = "#5BBAB5",
#                    "rickard-rickard_empirical" = "#589dd3")
COMMON.THEME <- theme(axis.title = element_text(size = 14),
                      axis.text = element_text(size = 12),
                      legend.title = element_text(size = 12),
                      legend.text = element_text(size = 12))

# Step 1. Load evaluation data
df.evaluation.metrics <- data.frame()
for (tsv.file in Sys.glob(paste0(EVALUATION.RESULTS.DIR, "/*evaluation_metrics.tsv"))) {
  df.temp <- read.csv(tsv.file, sep = "\t")
  df.temp$group <- paste0(df.temp$configuration_method, "-", df.temp$deconvolution_method)
  df.evaluation.metrics <- rbind(df.evaluation.metrics, df.temp)
}

# Step 2. Drop any group with fewer than 20 simulations
for (num.peptides in unique(df.evaluation.metrics$num_peptides)) {
  for (group in unique(df.evaluation.metrics$group)) {
    for (num.immunogenic.peptides in unique(df.evaluation.metrics$num_immunogenic_peptides)) {
      df.matched <- df.evaluation.metrics[
        (df.evaluation.metrics$num_peptides == num.peptides) &
        (df.evaluation.metrics$group == group) & 
        (df.evaluation.metrics$num_immunogenic_peptides == num.immunogenic.peptides),
      ]
      if (nrow(df.matched) < 20) {
        df.evaluation.metrics <- df.evaluation.metrics[
          !(df.evaluation.metrics$experiment_id %in% df.matched$experiment_id),
        ]
      }
    }
  }
}

# Step 3. Append DeconvoluteThis rows
df.deconvolutethis <- read.csv(DECONVOLUTETHIS.DESIGN.FILE)
df.deconvolutethis <- df.deconvolutethis[df.deconvolutethis$num_coverage == 3,]
for (num.peptides in unique(df.evaluation.metrics$num_peptides)) {
  for (num.immunogenic.peptides in unique(df.evaluation.metrics$num_immunogenic_peptides)) {
    df.matched <- df.deconvolutethis[
      (df.deconvolutethis$num_peptides == num.peptides) &
      (df.deconvolutethis$num_true_positive_peptides == num.immunogenic.peptides),
    ]
    num.peptides.per.pool <- df.matched$num_peptides_per_pool[1]
    num.total.pools <- df.matched$num_deconvolutethis_total_pools[1]
    num.init.pools <- (num.peptides / num.peptides.per.pool) * 3
    num.predicted.positive.peptides <- num.total.pools - num.init.pools
    precision <- num.immunogenic.peptides / num.predicted.positive.peptides
    sensitivity <- 1.0
    df.temp <- data.frame(
      experiment_id = c(paste0("deconvolutethis_", num.peptides,"peptides")),
      predicted_total_pools = c(num.total.pools),
      precision = c(precision),
      sensitivity = c(sensitivity),
      specificity = c(NA),
      aucroc = c(NA),
      num_peptides = c(num.peptides),
      num_peptides_per_pool = c(num.peptides.per.pool),
      num_coverage = c(3),
      num_immunogenic_peptides = c(num.immunogenic.peptides),
      rep_id = c("rep1"),
      configuration_method = c("deconvolutethis"),
      deconvolution_method = c("deconvolutethis"),
      group = c("deconvolutethis-deconvolutethis")
    )
    df.evaluation.metrics <- rbind(df.evaluation.metrics, df.temp)
  }
}

# Step 4. Define functions
PlotBarPlot <- function(df.evaluation.metrics, metric, num.peptides, y.axis.title) {
  df.plot <- df.evaluation.metrics[df.evaluation.metrics$num_peptides == num.peptides,]
  df.plot <- df.plot %>%
    dplyr::group_by(group, num_immunogenic_peptides) %>%
    dplyr::summarise(
      sd = sd(get(metric)),
      metric = mean(get(metric))
    )
  df.plot$num_immunogenic_peptides <- as.character(df.plot$num_immunogenic_peptides)
  for (group in unique(df.plot$group)) {
    for (num.immunogenic.peptides in c("1","2","3","4","5","6","7","8","9","10","15","20","25","30","40")) {
      df.matched <- df.plot[(df.plot$group == group) & (df.plot$num_immunogenic_peptides == num.immunogenic.peptides),]
      if (nrow(df.matched) == 0) {
        df.temp <- data.frame(
          group = c(group),
          num_immunogenic_peptides = c(num.immunogenic.peptides),
          sd = c(0.0),
          metric = c(0.0)
        )
        df.plot <- rbind(df.plot, df.temp)
      }
    }
  }
  df.plot$num_immunogenic_peptides <- factor(
    df.plot$num_immunogenic_peptides,
    levels = c("1","2","3","4","5",
               "6","7","8","9","10",
               "15","20","25","30","40")
  )
  # df.plot$group <- factor(
  #   df.plot$group,
  #   levels = c("ace-ace-em",
  #              "randomized-ace-em",
  #              "rickard-rickard-background-subtracted",
  #              "rickard-rickard_empirical",
  #              "rickard-rickard_bayesian",
  #              "deconvolutethis-deconvolutethis",
  #              "repeated-ace-em")
  # )
  p <- ggplot(df.plot, aes(x = num_immunogenic_peptides, y = metric, fill = group)) +
    geom_col(width = BOXPLOT.WIDTH, colour = "black", position = "dodge") +
    geom_errorbar(aes(ymin = metric - sd, ymax = ifelse(metric + sd > 1, 1, metric + sd)), width = 0.2, position = position_dodge(BOXPLOT.WIDTH)) +
    xlab("Number of Positive Peptides") + ylab(y.axis.title) + ggtitle("End-to-end ELISpot Experiment Benchmark") +
    scale_y_continuous(limits = c(0,1), breaks = c(0.00, 0.25, 0.50, 0.75, 1.00), expand = c(0,0)) +
    # scale_fill_manual(values = SOLVER.COLORS) +
    theme_pubr() +
    theme(plot.title = element_text(size = 12, hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          strip.text = element_text(size = 12),
          strip.background = element_rect(fill = "white"))
  return(p)
}

# Step 5. Generate plots
plot.120.peptides.precision <- PlotBarPlot(df.evaluation.metrics = df.evaluation.metrics, metric = "precision", num.peptides = 120, y.axis.title = "Precision")
plot.120.peptides.recall <- PlotBarPlot(df.evaluation.metrics = df.evaluation.metrics, metric = "sensitivity", num.peptides = 120, y.axis.title = "Recall")
print(plot.120.peptides.precision)
print(plot.120.peptides.recall)

plot.800.peptides.precision <- PlotBarPlot(df.evaluation.metrics = df.evaluation.metrics, metric = "precision", num.peptides = 800, y.axis.title = "Precision")
plot.800.peptides.recall <- PlotBarPlot(df.evaluation.metrics = df.evaluation.metrics, metric = "sensitivity", num.peptides = 800, y.axis.title = "Recall")
print(plot.800.peptides.precision)
print(plot.800.peptides.recall)
