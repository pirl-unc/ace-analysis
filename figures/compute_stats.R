library(dplyr)


# Figure 3 Statistics
# Panel B
TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/02_reduced_designs/reduced_designs_experiment_results_merged.tsv"
df <- read.csv(TSV.FILE, sep = "\t")
df.temp <- df %>%
    dplyr::filter(num_peptides == 100)
print("100/10/3x")
for (solver in unique(df.temp$solver)) {
  df.temp.2 <- df.temp[df.temp$solver == solver,]
  print(paste0(solver, " average precision: ", mean(df.temp.2$precision_empirical)))
  print(paste0(solver, " average number of total pools: ", mean(df.temp.2$predicted_total_pools)))
}
df.temp <- df %>%
  dplyr::filter(num_peptides == 800)
print("800/25/3x")
for (solver in unique(df.temp$solver)) {
  df.temp.2 <- df.temp[df.temp$solver == solver,]
  print(paste0(solver, " average precision: ", mean(df.temp.2$precision_empirical)))
  print(paste0(solver, " average number of total pools: ", mean(df.temp.2$predicted_total_pools)))
}
# Panel C
TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/03_deconvolution_methods/deconvolution_methods_experiment_results_merged.tsv"
df <- read.csv(TSV.FILE, sep = "\t")
df$group <- paste0(df$num_peptides, "/", df$num_peptides_per_pool, "/", df$num_coverage)
df$perc_positive_peptide_sequences <- floor((df$num_positive_peptide_sequences / df$num_peptides) * 100)
for (group in unique(df$group)) {
  df.temp <- df[(df$group == group) & (df$solver == "ace_golfy_clusteroff_noextrapools"),]
  print(group)
  print(paste0("Mean EM AUCROC: ", mean(df.temp$aucroc_score_em)))
  print(paste0("Mean LASSO AUCROC: ", mean(df.temp$aucroc_score_lasso)))
  print(paste0("Mean Empirical AUCROC: ", mean(df.temp$aucroc_score_empirical)))
}
for (group in unique(df$group)) {
  df.temp <- df[
    (df$group == group) & 
    (df$perc_positive_peptide_sequences <= 3) &
    (df$solver == "ace_golfy_clusteroff_noextrapools"),
  ]
  print(group)
  print(paste0("Mean EM AUCROC (<= 5%) : ", mean(df.temp$aucroc_score_em)))
  print(paste0("Mean LASSO AUCROC (<= 5%): ", mean(df.temp$aucroc_score_lasso)))
  print(paste0("Mean Empirical AUCROC (<= 5%): ", mean(df.temp$aucroc_score_empirical)))
}


# Figure 4 Statistics
# Panel C and D
TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/05_alanine_scanning/held_out_data/alanine_scanning_experiment_results_merged.tsv"
df <- read.csv(TSV.FILE, sep = "\t")
for (num.peptides in unique(df$num_peptides)) {
  df.temp <- df[df$num_peptides == num.peptides,]
  print(paste0(num.peptides, " peptides"))
  for (solver in unique(df.temp$solver)) {
    df.temp.2 <- df.temp[df.temp$solver == solver,]
    print(paste0(solver, " mean precision = ", mean(df.temp.2$precision_empirical)))
  }
  print("")
}
# Panel E
TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/06_real_world_datasets/tan_et_al_inf_gen_evol_2021/noextrapools/177peptides_10perpool_3x/tan_et_al_inf_gen_evol_2021_experiment_results_585570736.tsv"
df <- read.csv(TSV.FILE, sep = "\t")
for (solver in unique(df$solver)) {
  df.temp <- df[df$solver == solver, ]
  print(paste0(solver, " mean precision = ", mean(df.temp$precision_empirical)))
  print(paste0(solver, " mean num total pools = ", mean(df.temp$predicted_total_pools)))
}
# Panel F
TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/06_real_world_datasets/cameron_et_al_sci_trans_med_2013/36peptides_6perpool_3x/cameron_et_al_sci_trans_med_2013_experiment_results_207124847.tsv"
df <- read.csv(TSV.FILE, sep = "\t")
for (solver in unique(df$solver)) {
  df.temp <- df[df$solver == solver, ]
  print(paste0(solver, " mean precision = ", mean(df.temp$precision_empirical)))
  print(paste0(solver, " mean num total pools = ", mean(df.temp$predicted_total_pools)))
}

# Supplementary Figure 5 Statistics
TSV.FILE <- "/Users/leework/Documents/Research/projects/project_ace/data/processed/01_benchmark_ace/05_alanine_scanning/training_data/alanine_scanning_experiment_results_merged.tsv"
df <- read.csv(TSV.FILE, sep = "\t")
for (num.peptides in unique(df$num_peptides)) {
  df.temp <- df[df$num_peptides == num.peptides,]
  print(paste0(num.peptides, " peptides"))
  for (solver in unique(df.temp$solver)) {
    df.temp.2 <- df.temp[df.temp$solver == solver,]
    print(paste0(solver, " mean precision = ", mean(df.temp.2$precision_empirical)))
  }
  print("")
}
