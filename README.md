
<!-- README.md is generated from README.Rmd. Please edit README.Rmd file -->

# `MiXcan: Statistical Framework for Cell-type-Aware Transcriptome-Wide Association Studies with Bulk Tissue Data`

## Introduction to **MiXcan**

**Goal:**

- Constructs cell-type-level prediction models for genetically regulated
  expression (GReX);

- Predicts cell-type-level GReX in new genotype data; and

- Performs cell-type-aware TWAS.

**Advantages over tissue-level TWAS:**

- Improves GReX prediction accuracy;

- Boosts the study power, especially for genes that function in minor
  cell types or have different association directions in different cell
  types;

- Sheds light on the responsible cell type(s) of associations.

**Disadvantages over tissue-level TWAS:**

- Requires prior knowledge on disease-critical cell types and their
  proportions in tissue;

- Has more model parameters;

- May be less powerful than tissue-level TWAS for genes that have
  similar disease associations in different cell types or function in
  major cell types.

**Input:**

- Prediction model construction: genotype, covariates, and gene
  expression data (same as in PrediXcan) + cell-type composition
  estimates (e.g. from existing methods, such as ESTIMATE, CIBERSORT,
  xCell).

- Association Analysis: genotype, covariates and phenotype data (same as
  in PrediXcan).

**Output:**

- Prediction model construction: Cell-type-specific or nonspecific
  prediction weights for different genes.

- Association Analysis: Tissue-level association p-values and
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
devtools::install_github("songxiaoyu/MiXcan/Package")
```

The typical install time of the package is less than 5 minutes.

## Example of use

Below demonstrates the MiXcan analysis pipeline on a single peusdo gene.
In reality, multiple genes can be analyzed in parallel. Holdout set can
be pre-excluded to allow model training on the remaining samples.

### Data

The peusdo data are included in the Github page. We can load the data:

``` r
library(MiXcan)
data(example_data)
```

### MiXcan analysis pipeline

Step 1: Improving the estimation of the cell-type composition from
prior. This step is optional, and can be ignored if input cell-type
composition estimates is preferred.

``` r
set.seed(123)
pi_estimation_result <- pi_estimation(expression_matrix = GTEx_epithelial_genes,
              prior = GTEx_prior, 
              n_iteration = 5) 
pi_example=pi_estimation_result$mean_trim_0.05
pi_example[1:3]
```

    ## [1] 0.2394073 0.2258077 0.1000000

``` r
# Note: 0<prior<1 that prior should not hit the boundaries. If xCell score is used, add a small perturbation for zero.
```

Step 2: Estimating cell-type-specific (and nonspecific) GReX prediction
weights of a gene using the MiXcan function

``` r
set.seed(111)
foldid_example <- sample(1:10, length(y_example), replace=T)
MiXcan_result <- MiXcan(y=y_example, x=x_example, cov = cov_example, pi= pi_example, foldid = foldid_example)
MiXcan_result$beta.SNP.cell1
```

    ##    xNameMatrix      weight
    ## 1         SNP1  0.00000000
    ## 2         SNP2  0.00000000
    ## 3         SNP3  0.01675598
    ## 4         SNP4  0.05282617
    ## 5         SNP5  0.00000000
    ## 6         SNP6  0.00000000
    ## 7         SNP7  0.00000000
    ## 8         SNP8  0.00000000
    ## 9         SNP9  0.00000000
    ## 10       SNP10 -0.01754009
    ## 11       SNP11  0.00000000
    ## 12       SNP12  0.00000000
    ## 13       SNP13  0.00000000
    ## 14       SNP14  0.00000000
    ## 15       SNP15  0.00000000
    ## 16       SNP16  0.00000000
    ## 17       SNP17  0.00000000
    ## 18       SNP18  0.00000000
    ## 19       SNP19 -0.07213873

``` r
MiXcan_result$beta.SNP.cell2
```

    ##    xNameMatrix      weight
    ## 1         SNP1  0.00000000
    ## 2         SNP2  0.00000000
    ## 3         SNP3  0.01675598
    ## 4         SNP4  0.05282617
    ## 5         SNP5  0.00000000
    ## 6         SNP6  0.00000000
    ## 7         SNP7  0.00000000
    ## 8         SNP8  0.00000000
    ## 9         SNP9  0.00000000
    ## 10       SNP10 -0.01754009
    ## 11       SNP11  0.00000000
    ## 12       SNP12  0.00000000
    ## 13       SNP13  0.00000000
    ## 14       SNP14  0.00000000
    ## 15       SNP15  0.00000000
    ## 16       SNP16  0.00000000
    ## 17       SNP17  0.00000000
    ## 18       SNP18  0.00000000
    ## 19       SNP19 -0.07213873

Step 3: Extracting the weights from the MiXcan output.

``` r
MiXcan_weight_result <- MiXcan_extract_weight(model = MiXcan_result)
```

    ## Joining with `by = join_by(xNameMatrix)`

``` r
MiXcan_weight_result
```

    ##   xNameMatrix weight_cell_1 weight_cell_2        type
    ## 1        SNP3    0.01675598    0.01675598 NonSpecific
    ## 2        SNP4    0.05282617    0.05282617 NonSpecific
    ## 3       SNP10   -0.01754009   -0.01754009 NonSpecific
    ## 4       SNP19   -0.07213873   -0.07213873 NonSpecific

Step 4: Predicting the cell-type-specific or nonspecific expression
levels of the gene in a new genetic data.

``` r
MiXcan_prediction_result <- MiXcan_prediction(weight = MiXcan_weight_result, new_x = new_X_example)
MiXcan_prediction_result[1:3,]
```

    ##          cell_1     cell_2
    ## [1,] 0.00000000 0.00000000
    ## [2,] 0.03450198 0.03450198
    ## [3,] 0.05204207 0.05204207

Step 5: Association analysis with MiXcan predicted GReX levels

``` r
MiXcan_association_result <- MiXcan_association(new_y = MiXcan_prediction_result,
                                                new_cov = new_cov_example, 
                                                new_outcome = new_outcome_example, family  = "binomial")
MiXcan_association_result
```

    ##    cell1_est cell1_se   cell1_p  cell2_est cell2_se   cell2_p p_combined
    ## 1 -0.8324663 1.697274 0.6237995 -0.8324663 1.697274 0.6237995  0.6237995

## Pretrained models:

Pretrained models are provided here for epithelial vs. stromal cell
types in mammary tissues. A total of 125 European ancestry women in GTEx
was used for model trainig. Training data can be accessed from dbGap
(Study Accession: phs000424.v9.p2).

Note: As the proof of concept, a subset of SNPs reported in the mammary
tissues in predictdb (<https://predictdb.org>) were used as the SNP pool
for model training.

File Path:
“PreTrainedModel/MiXcan_model_weights_trained_in_GTEx_v8_mammary.tsv”

Users of the pre-trained models can apply the weights to new genotype
data as in Step 4-5 for cell-type-aware transcriptome-wide association
studies.
