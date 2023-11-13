library(dplyr)
library(stringr)


OUTPUT.DIR <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/tables"


# 01. Main Table 2
EVALUATION.RESULTS.DIR <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/02_benchmark_ace/01_benchmark_against_existing_tools/06_evaluation_results/300immunogenic_30nonimmunogenic_1dispersion_fnr0.00"
df.evaluation.metrics <- data.frame()
for (tsv.file in Sys.glob(paste0(EVALUATION.RESULTS.DIR, "/*evaluation_metrics.tsv"))) {
  df.temp <- read.csv(tsv.file, sep = "\t")
  df.temp$group <- paste0(df.temp$configuration_method, "_", df.temp$deconvolution_method)
  df.evaluation.metrics <- rbind(df.evaluation.metrics, df.temp)
}

df.precision.recall.aucroc.numpools <- df.evaluation.metrics %>%
  dplyr::filter(configuration_method %in% c("ace-s","ace","randomized","repeated")) %>%
  dplyr::filter(num_peptides == 120) %>%
  dplyr::group_by(configuration_method, deconvolution_method) %>%
  dplyr::summarise(aucroc = mean(aucroc),
                   precision = mean(precision),
                   sensitivity = mean(sensitivity),
                   num_total_pools = mean(predicted_total_pools))

df.evaluation.results <- data.frame()
for (tsv.file in Sys.glob(paste0(EVALUATION.RESULTS.DIR, "/*evaluation_results.tsv"))) {
  df.temp <- read.csv(tsv.file, sep = "\t")
  df.evaluation.results <- rbind(df.evaluation.results, df.temp)
}
df.evaluation.results <- df.evaluation.results %>%  
  dplyr::filter(configuration_method %in% c("ace-s","ace","randomized","repeated")) %>%
  dplyr::filter(num_peptides == 120)
df.evaluation.results$group <- paste0(df.evaluation.results$configuration_method, "_", df.evaluation.results$deconvolution_method)
df.aucprc <- data.frame()
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
        (df.evaluation.results$num_immunogenic_peptides == num.immunogenic.peptides),
    ]
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
  df.auc.ci.prc$num_immunogenic_peptides <- num.immunogenic.peptides
  df.aucprc <- rbind(df.aucprc, df.auc.ci.prc)
}

df.aucprc.all <- df.aucprc %>%
  dplyr::group_by(modnames) %>%
  dplyr::summarise(aucprc = mean(mean))
split_values <- str_split(df.aucprc.all$modnames, "_", simplify = TRUE)
df.aucprc.all$configuration_method <- split_values[, 1]
df.aucprc.all$deconvolution_method <- split_values[, 2]

df.main.table.2 <- merge(df.precision.recall.aucroc.numpools, df.aucprc.all, by = c("configuration_method", "deconvolution_method"))
write.table(x = df.main.table.2, file = paste0(OUTPUT.DIR, "/main_table_2.tsv"), sep = "\t", row.names = F)
