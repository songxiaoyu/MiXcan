
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `MiXcan: Statistical Framework for Cell-type-Specific Transcriptome-Wide Association Studies with Bulk Tissue Data`

## Introduction to **MiXcan**

**Goal:**

-   Construct cell-type-specific gene expression prediction models;

-   Apply these models to predict cell-type-specific gene expression
    levels in new genotype data; and

-   Perform cell-type-specific TWAS.

**Advantages over tissue-level TWAS:**

-   Improve expression prediction accuracy;

-   Boost the study power, especially for genes that function in minor
    cell types or have different association directions in different
    cell types;

-   Shed light on the responsible cell type(s) of associations.

**Disadvantages over tissue-level TWAS:**

-   Require prior knowledge on cell types

-   Increased complexity with more model parameters

-   May be less powerful for genes that (1) have different associations
    with genotypes in different cell types, and (2) have similar
    associations with phenotypes in different cell types or are
    associated with disease in major cell types.

**Input:**

-   Prediction model construction: Same as in PrediXcan (genotype,
    covariates, and gene expression data) + prior cell-type composition
    estimates (e.g. from existing methods, such as ESTIMATE, CIBERSORT,
    xCell).

-   Association Analysis: Same as in PrediXcan (genotype, covariates and
    phenotype data).

**Output:**

-   Prediction model construction: Cell-type-specific and nonspecific
    prediction weights.

-   Association Analysis: Tissue-level association p-values and
    cell-type-level association summaries including estimates, standard
    error and p-values.

A full description of the method can be found in our
[paper](https://www.biorxiv.org/content/10.1101/2022.03.15.484509v1.abstract).

``` r
knitr::opts_chunk$set(echo = TRUE)
library(MiXcan) 
```

## Installation

### Hardware Requirements

The MiXcan package requires only a standard computer with a reasonable
RAM/CPU to support the operations. The minimal RAM recommended by
authors is 2 GB. Authors used a computer with the following
specification:

RAM: 32 GB

CPU: 2.3 GHz 8-core

### Software Requirements

The github package is developed on macOS operating systems. The MiXcan
pacakge is developed under an open source software R (version 4.1.2).
Different versions of R can be downloaded at
<https://cran.r-project.org/>.

### Package Installation

With R, users can install the MiXcan package directly from GitHub with
[devtools](https://github.com/hadley/devtools):

``` r
install.packages("devtools")
devtools::install_github("songxiaoyu/MiXcan")
```

The typical install time of the package is less than 5 minutes.

## Example of use

Below demonstrates the MiXcan analysis pipeline on a single peusdo gene.
In reality, multiple genes can be analyzed in parallel.

### Data

The peusdo data are included in the Github page. We can load the data:

``` r
library(MiXcan)
data(example_data)
```

### MiXcan analysis pipeline

Step 1 (optional): Improving the estimation of the cell-type composition
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

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
    ## ✓ tibble  3.1.6     ✓ dplyr   1.0.8
    ## ✓ tidyr   1.2.0     ✓ stringr 1.4.0
    ## ✓ readr   2.1.2     ✓ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x purrr::accumulate() masks foreach::accumulate()
    ## x dplyr::filter()     masks stats::filter()
    ## x dplyr::lag()        masks stats::lag()
    ## x purrr::when()       masks foreach::when()

``` r
nCores=detectCores()-1; registerDoParallel(nCores) # use parallel computing but leave 1 core out. 
# Note: 0<prior<1 that prior should not hit the boundaries. If xCell score is used, check zero.
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

Step 2: Estimating cell-type-specific (and nonspecific) prediction
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

Step 3: Extract the weights from the output of MiXcan function.

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

Step 4: Predict the cell-type-specific or nonspecific expression levels
of a gene with MiXcan model in new genetic data.

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

Step 5: Association analysis with MiXcan predicted gene expression
levels

``` r
MiXcan_association_result <- MiXcan_association(MiXcan_predicted_expr = MiXcan_prediction_result,
                                                covariates = covariates_example, outcome = outcome_example, family  = "binomial")
MiXcan_association_result
```

    ##   cell_1_estimate cell_1_std.error  cell_1_p cell_2_estimate cell_2_std.error
    ## 1      -0.1465637         1.778105 0.9343073      -0.1465637         1.778105
    ##    cell_2_p p_combined
    ## 1 0.9343073  0.9343073

## Pretrained models in mammary tissues

Data for training the cell-type-specific (and nonspecific) prediction
models are protected and can be accessed via dbGap. Pre-trained
prediction models (Step 1-3) using the mammary tissues from 125 Eurpean
Ancestry (EA) can be accessed by

``` r
weights=read.table("data/MiXcan_model_weights_trained_in_GTEx_v8_mammary.tsv", header=T)
```

Users of the pre-trained models can apply the weights to new genotype
data as in Step 4-5 for cell-type-specific association analyses.
