library(precrec)


# 01. Define constants
EVALUATION.RESULTS.DIR <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/02_benchmark_ace/01_benchmark_against_existing_tools/06_evaluation_results/300immunogenic_5nonimmunogenic_1dispersion_fnr0.00"
OUTPUT.DIR <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/02_benchmark_ace/01_benchmark_against_existing_tools/06_evaluation_results/300immunogenic_5nonimmunogenic_1dispersion_fnr0.00"

# 02. Load data
df.evaluation.results <- data.frame()
for (tsv.file in Sys.glob(paste0(EVALUATION.RESULTS.DIR, "/*evaluation_results.tsv"))) {
  df.temp <- read.csv(tsv.file, sep = "\t")
  df.temp$group <- paste0(df.temp$configuration_method, "_", df.temp$deconvolution_method)
  df.evaluation.results <- rbind(df.evaluation.results, df.temp)
}

# 03. Compute AUC
df.auc.ci.all <- data.frame()
for (num.peptides in c(120, 800)) {
  for (num.immunogenic.peptides in c(1,2,3,4,5,6,7,8,9,10,15,20,25,30,40)) {
    list.scores <- list()
    list.labels <- list()
    dsids <- c()
    mod.names <- c()
    idx <- 1
    dsid <- 1
    for (group in unique(df.evaluation.results$group)) {
      df.temp <- df.evaluation.results[
        (df.evaluation.results$num_peptides == num.peptides) &
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
    df.auc.ci$numpeptides <- num.peptides
    df.auc.ci$numimmunogenicpeptides <- num.immunogenic.peptides
    df.auc.ci.all <- rbind(df.auc.ci.all, df.auc.ci)
  }
}

write.table(x = df.auc.ci.all, 
            file = paste0(OUTPUT.DIR, "/ablation_study.tsv"), 
            sep = "\t",
            row.names = F)
