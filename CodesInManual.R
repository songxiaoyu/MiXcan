rm(list=ls())
library(MiXcan)
library(doParallel)
library(tidyverse)
nCores=detectCores()-1; registerDoParallel(nCores) # use parallel computing for speed,

save(GTEx_epithelial_genes, GTEx_prior,
     x_example, y_example, cov_example,
     new_X_example, new_outcome_example, new_cov_example,
     file="example_data.rda")
a=colnames(new_X_example)
new_X_example=matrix(rbinom(2000, 2, 0.3), ncol=4)
colnames(new_X_example)=c("SNP3", "SNP4", "SNP10", "SNP19")
new_outcome_example=rbinom(500, 1, 0.3)
new_cov_example=cbind(age=round(rnorm(500, 50), 2), female=rbinom(500, 1, 0.45))


data(example_data)
# pi est
set.seed(123)
pi_estimation_result <- pi_estimation(expression_matrix = GTEx_epithelial_genes,
                                      prior = GTEx_prior,
                                      n_iteration = 5)
pi_example=pi_estimation_result$mean_trim_0.05
pi_example[1:3]

# MiXcan fit
set.seed(111)
foldid_example <- sample(1:10, length(y_example), replace=T)

MiXcan_result <- MiXcan(y=y_example, x=x_example, cov = cov_example,
                        pi= pi_example,
                        foldid = foldid_example, yName="Gene1")
MiXcan_result$beta.SNP.cell1
MiXcan_result$beta.SNP.cell2

# Extract Info
MiXcan_weight_result <- MiXcan_extract_weight(model = MiXcan_result)
MiXcan_weight_result

MiXcan_summary_result <- MiXcan_extract_summary(x=x_example, y=y_example,
                                                pi=pi_example, model=MiXcan_result)
MiXcan_summary_result

# Refit
MiXcan_weight_refit <- MiXcan_refit_weight(model = MiXcan_result,
                                           y=y_example,
                                           x=x_example, cov = cov_example,
                                           pi= pi_example)
MiXcan_weight_refit
MiXcan_summary_refit <- MiXcan_refit_summary(model = MiXcan_result, y=y_example,
                                             x=x_example, cov = cov_example,
                                             pi= pi_example)

MiXcan_summary_refit

# Prediction in a new data
MiXcan_prediction_result <- MiXcan_prediction(weight = MiXcan_weight_result, new_x = new_X_example)
MiXcan_prediction_result

# Association in a new data
MiXcan_association_result <- MiXcan_association(MiXcan_predicted_expr = MiXcan_prediction_result,
                                                covariates = covariates_example, outcome = outcome_example, family  = "binomial")
MiXcan_association_result
