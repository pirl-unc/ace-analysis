library(dplyr)
library(stringr)


OUTPUT.DIR <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/tables"


# 01. Main Table 2
EVALUATION.RESULTS.DIR <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/02_benchmark_ace/01_benchmark_against_existing_tools/06_evaluation_results/300immunogenic_30nonimmunogenic_1dispersion_fnr0.00"
DECONVOLUTETHIS.DESIGN.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/raw/references/deconvolutethis_designs.csv"
MIN.SIMULATIONS <- 20
df.evaluation.metrics <- data.frame()
for (tsv.file in Sys.glob(paste0(EVALUATION.RESULTS.DIR, "/*evaluation_metrics.tsv"))) {
  df.temp <- read.csv(tsv.file, sep = "\t")
  df.temp$group <- paste0(df.temp$configuration_method, "_", df.temp$deconvolution_method)
  df.evaluation.metrics <- rbind(df.evaluation.metrics, df.temp)
}

# Append DeconvoluteThis rows
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
      group = c("deconvolutethis")
    )
    df.evaluation.metrics <- rbind(df.evaluation.metrics, df.temp)
  }
}
# Step 2. Drop any group with fewer than 20 simulations
for (num.peptides in unique(df.evaluation.metrics$num_peptides)) {
  for (group in unique(df.evaluation.metrics$group)) {
    if (group == "deconvolutethis"){
      next
    }
    for (num.immunogenic.peptides in unique(df.evaluation.metrics$num_immunogenic_peptides)) {
      df.matched <- df.evaluation.metrics[
        (df.evaluation.metrics$num_peptides == num.peptides) &
          (df.evaluation.metrics$group == group) & 
          (df.evaluation.metrics$num_immunogenic_peptides == num.immunogenic.peptides),
      ]
      if (nrow(df.matched) < MIN.SIMULATIONS) {
        df.evaluation.metrics <- df.evaluation.metrics[
          !(df.evaluation.metrics$experiment_id %in% df.matched$experiment_id),
        ]
      }
    }
  }
}
df.precision.recall.aucroc.numpools <- df.evaluation.metrics %>%
  dplyr::filter(group %in% c("ace-s_ace-cem", 
                             "randomized_ace-cem",
                             "repeated_ace-empirical",
                             "strandberg_strandberg-background-subtracted",
                             "strandberg_strandberg_bayesian",
                             "strandberg_strandberg-empirical",
                             "deconvolutethis")) %>%
  dplyr::group_by(configuration_method, deconvolution_method, num_peptides) %>%
  dplyr::summarise(aucroc = mean(aucroc),
                   precision = mean(precision),
                   sensitivity = mean(sensitivity),
                   num_total_pools = mean(predicted_total_pools))

df.evaluation.results <- data.frame()
for (tsv.file in Sys.glob(paste0(EVALUATION.RESULTS.DIR, "/*evaluation_results.tsv"))) {
  df.temp <- read.csv(tsv.file, sep = "\t")
  df.temp$group <- paste0(df.temp$configuration_method, "_", df.temp$deconvolution_method)
  df.evaluation.results <- rbind(df.evaluation.results, df.temp)
}
df.evaluation.results <- df.evaluation.results %>%
  dplyr::filter(group %in% c("ace-s_ace-cem",
                             "randomized_ace-cem",
                             "repeated_ace-empirical",
                             "strandberg_strandberg-background-subtracted",
                             "strandberg_strandberg_bayesian",
                             "strandberg_strandberg-empirical"))
df.aucprc <- data.frame()
for (num.peptides in unique(df.evaluation.results$num_peptides)) {
  for (num.immunogenic.peptides in unique(df.evaluation.results$num_immunogenic_peptides)) {
    list.scores <- list()
    list.labels <- list()
    dsids <- c()
    mod.names <- c()
    idx <- 1
    dsid <- 1
    for (group in unique(df.evaluation.results$group)) {
      df.temp <- df.evaluation.results[
        (df.evaluation.results$group == group) &
        (df.evaluation.results$num_peptides == num.peptides) &
        (df.evaluation.results$num_immunogenic_peptides == num.immunogenic.peptides),
      ]
      if (length(unique(df.temp$rep_id)) < MIN.SIMULATIONS) {
        next
      }
      for (rep.id in unique(df.temp$rep_id)) {
        df.temp_ <- df.temp[df.temp$rep_id == rep.id,]
        list.scores[[idx]] <- df.temp_$peptide_spot_count
        list.labels[[idx]] <- df.temp_$binding
        dsids <- c(dsids, dsid)
        mod.names <- c(mod.names, group)
        idx <- idx + 1
      }
      dsid <- dsid + 1
    }
    smmdat <- mmdata(scores = list.scores, 
                     labels = list.labels, 
                     modnames = mod.names, 
                     dsids = dsids)
    smcurves <- evalmod(smmdat, cb_alpha = 0.05)
    df.auc.ci <- auc_ci(smcurves)
    df.auc.ci.prc <- df.auc.ci[df.auc.ci$curvetypes == "PRC",]
    df.auc.ci.prc$num_peptides <- num.peptides
    df.auc.ci.prc$num_immunogenic_peptides <- num.immunogenic.peptides
    df.aucprc <- rbind(df.aucprc, df.auc.ci.prc)
  }
}

df.aucprc.all <- df.aucprc %>%
  dplyr::group_by(modnames) %>%
  dplyr::summarise(aucprc = mean(mean))
split_values <- str_split(df.aucprc.all$modnames, "_", simplify = TRUE)
df.aucprc.all$configuration_method <- split_values[, 1]
df.aucprc.all$deconvolution_method <- split_values[, 2]

df.main.table.2 <- merge(df.precision.recall.aucroc.numpools, df.aucprc.all, by = c("configuration_method", "deconvolution_method"))
write.table(x = df.main.table.2, file = paste0(OUTPUT.DIR, "/main_table_2.tsv"), sep = "\t", row.names = F)
