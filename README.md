
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `MiXcan: Statistical Models for Cell-type Specific Transcriptome-Wide Association Studies with Bulk Tissue Data`

## Introduction

Considering the cell-type composition of a tissue, the goal of
**MiXcan** is to

  - Provide the cell-type specific gene expression levels and improve
    the prediction accuracy for the tissue.

  - Boost the study power for gene identifications in TWAS and shed
    light on the functional cell type(s) of the associations.

A full description of the method can be found in our
[paper](url%20link).

``` r
knitr::opts_chunk$set(echo = TRUE)
library(MiXcan)
```

### Installation

You can install the latest version directly from GitHub with
[devtools](https://github.com/hadley/devtools):

``` r
install.packages("devtools")
devtools::install_github("songxiaoyu/MiXcan")
```

## Example of use

Below is an example of MiXcan analysis pipeline.

### Data

The sample data are included in the Github page. We will load the data:

``` r
library(MiXcan)
load("data/example_data.rda")
```

### MiXcan analysis pipeline

Step 1 (option): Improving the estimation of the cell-type composition
Pi.

``` r
library(doParallel)
```

    ## Loading required package: foreach

    ## Loading required package: iterators

    ## Loading required package: parallel

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.3     ✓ purrr   0.3.4
    ## ✓ tibble  3.0.4     ✓ dplyr   1.0.7
    ## ✓ tidyr   1.1.2     ✓ stringr 1.4.0
    ## ✓ readr   1.4.0     ✓ forcats 0.5.0

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x purrr::accumulate() masks foreach::accumulate()
    ## x dplyr::filter()     masks stats::filter()
    ## x dplyr::lag()        masks stats::lag()
    ## x purrr::when()       masks foreach::when()

``` r
nCores=detectCores()-1; registerDoParallel(nCores) # use parallel computing but leave 1 core out. 
pi_estimation_result <- pi_estimation(expression_matrix = GTEx_epithelial_genes,
              prior = GTEx_prior,
              n_iteration = 5)
```

    ## [1] "EM algorithm converged"
    ## [1] 1
    ## [1] "EM algorithm converged"
    ## [1] 2
    ## [1] "EM algorithm converged"
    ## [1] 3
    ## [1] "EM algorithm converged"
    ## [1] 4
    ## [1] "EM algorithm converged"
    ## [1] 5

Step 2. Estimating cell-type specific (and non-specific) prediction
weights for the expression levels of a gene using the MiXcan function

``` r
set.seed(111)
foldid_example <- sample(1:10, length(y_example), replace=T)
MiXcan_result <- MiXcan(y=y_example, x=x_example, cov = cov_example, pi= pi_estimation_result$mean_trim_0.05, foldid = foldid_example)
```

    ## [1] 39
    ##       2.5%      97.5% 
    ## -0.8464646  2.3792982 
    ## [1] "NonSpecific"

``` r
MiXcan_result$beta.SNP.cell1
```

    ##    nameMatrix      weight
    ## 1        SNP1  0.00000000
    ## 2        SNP2  0.00000000
    ## 3        SNP3  0.00614620
    ## 4        SNP4  0.04108504
    ## 5        SNP5  0.00000000
    ## 6        SNP6  0.00000000
    ## 7        SNP7  0.00000000
    ## 8        SNP8  0.00000000
    ## 9        SNP9  0.00000000
    ## 10      SNP10 -0.03021309
    ## 11      SNP11  0.00000000
    ## 12      SNP12  0.00000000
    ## 13      SNP13  0.00000000
    ## 14      SNP14  0.00000000
    ## 15      SNP15  0.00000000
    ## 16      SNP16  0.00000000
    ## 17      SNP17  0.00000000
    ## 18      SNP18  0.00000000
    ## 19      SNP19 -0.06553891

``` r
MiXcan_result$beta.SNP.cell2
```

    ##    nameMatrix      weight
    ## 1        SNP1  0.00000000
    ## 2        SNP2  0.00000000
    ## 3        SNP3  0.00614620
    ## 4        SNP4  0.04108504
    ## 5        SNP5  0.00000000
    ## 6        SNP6  0.00000000
    ## 7        SNP7  0.00000000
    ## 8        SNP8  0.00000000
    ## 9        SNP9  0.00000000
    ## 10      SNP10 -0.03021309
    ## 11      SNP11  0.00000000
    ## 12      SNP12  0.00000000
    ## 13      SNP13  0.00000000
    ## 14      SNP14  0.00000000
    ## 15      SNP15  0.00000000
    ## 16      SNP16  0.00000000
    ## 17      SNP17  0.00000000
    ## 18      SNP18  0.00000000
    ## 19      SNP19 -0.06553891

3.  Extract the weights from the output of MiXcan
function.

<!-- end list -->

``` r
MiXcan_weight_result <- MiXcan_extract_weight(MiXcan_model = MiXcan_result)
```

    ## Joining, by = "nameMatrix"

``` r
MiXcan_weight_result
```

    ##   nameMatrix weight_cell_1 weight_cell_2        type
    ## 1       SNP3    0.00614620    0.00614620 NonSpecific
    ## 2       SNP4    0.04108504    0.04108504 NonSpecific
    ## 3      SNP10   -0.03021309   -0.03021309 NonSpecific
    ## 4      SNP19   -0.06553891   -0.06553891 NonSpecific

4.  Predict the cell-type specific or non-specific expression levels of
    a gene with MiXcan model in new genetic
data.

<!-- end list -->

``` r
MiXcan_prediction_result <- MiXcan_prediction(weight = MiXcan_weight_result, new_x = new_X_example)
MiXcan_prediction_result
```

    ##           cell_1      cell_2
    ## id1  -0.05324651 -0.05324651
    ## id2   0.00000000  0.00000000
    ## id3  -0.03150262 -0.03150262
    ## id4  -0.05466697 -0.05466697
    ## id5  -0.12020588 -0.12020588
    ## id6  -0.08345960 -0.08345960
    ## id7  -0.12596509 -0.12596509
    ## id8  -0.17921161 -0.17921161
    ## id9  -0.07770039 -0.07770039
    ## id10 -0.02406689 -0.02406689

5.  Association analysis with MiXcan predicted gene expression
levels

<!-- end list -->

``` r
MiXcan_association_result <- MiXcan_association(MiXcan_predicted_expr = MiXcan_prediction_result,
                                                covariates = covariates_example, outcome = outcome_example, family  = "binomial")
MiXcan_association_result
```

    ##   cell_1_estimate cell_1_std.error  cell_1_p cell_2_estimate cell_2_std.error
    ## 1      -0.1465637         1.778105 0.9343073      -0.1465637         1.778105
    ##    cell_2_p p_combined
    ## 1 0.9343073  0.9343073
