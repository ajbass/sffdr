
# sffdr package

<!-- badges: start -->

[![R-CMD-check](https://github.com/ajbass/sffdr/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/ajbass/sffdr/actions/workflows/R-CMD-check.yml)

<!---[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/sffdr)](https://cran.r-project.org/package=lit)--->
<!-- badges: end -->

## Overview

<img src="inst/figures/sffdr.png" align="right" height="150" />

The `sffdr` package implements the surrogate functional false discovery
rate (sfFDR) procedure which integrates GWAS summary statistics of
related traits to increase power within the fFDR framework. The inputs
into `sffdr` are a set of p-values from a GWAS of interest and a set of
p-values (or test statistics) from one or many informative GWAS. There
are a few key quantities estimated by the package: the functional
q-value, functional local FDR, and functional p-value.

The quantities estimated by `sffdr` can be used for a variety of
significance analyses in genome-wide studies while incorporating
informative data: the functional q-values provide positive FDR (related
to FDR) control, the functional local FDRs can be used for post-GWAS
analysis such as functional fine mapping, and the functional p-values
provide type I error rate control.

See our manuscript for additional details:

> Bass AJ, Wallace C. Exploiting pleiotropy to enhance variant discovery
> with functional false discovery rates.
> [medRxiv](https://www.medrxiv.org/content/10.1101/2024.09.24.24314276v1);
> 2024.

Note that this work is an extension of the functional FDR methodology
and the software builds on some of the functions in the original `fFDR`
package found at <https://github.com/StoreyLab/fFDR>.

## Installation

This software is implemented in the `R` statistical programming
language. The development version of `sffdr` can be installed using the
following code:

``` r
# install devtools
install.packages("devtools")
devtools::install_github("ajbass/sffdr")
```

The vignette can be viewed by typing:

``` r
browseVignettes(package = "sffdr")
```

## Quick start

Load the `sffdr` package and example data set:

``` r
library(sffdr)
set.seed(123)
data(bmi)
head(sumstats)
p <- sumstats$bmi
z <- as.matrix(sumstats[,-1])
```

The `sumstats` is a data frame contains 10,000 independent p-values for
body mass index (BMI), body fat percentage (BFP), cholesterol, and
triglycerides. Our primary trait of interest is BMI and the informative
traits are BFP, cholesterol and triglycerides.

We first model the relationship between the functional proportion of
true null hypotheses and the summary statistics using the `fpi0est`
function. The function `pi0_model` can be used to help create the design
matrix for `fpi0est`:

``` r
# Create model: choose knots at small quantiles
mpi0 <- pi0_model(z = z, knots = c(0.01, 0.025, 0.05, 0.1))

# Estimation of functional pi0 using the design matrix in mpi0
fpi0 <- fpi0est(p = p,
                z = mpi0$zt,
                pi0_model = mpi0$fmod)
```

The functional FDR and p-value quantities can be estimated with the
`sffdr` function:

``` r
# apply sfFDR
sffdr_out <- sffdr(p, fpi0 = fpi0$fpi0, epsilon = min(p))   

# plot significance results
plot(sffdr_out)

# Functional P-values, Q-values, and local FDR
fp <- sffdr_out$fpvalues
fq <- sffdr_out$fqvalues
flfdr <- sffdr_out$flfdr
```

The functional q-value (`fqvalue`) is a measure of significance in terms
of the positive FDR (closely related to FDR), the local FDR (`flfdr`) is
a posterior error probability, and the functional p-value (`fpvalue`)
can be interpreted similarly to a standard p-value.

Thus far, we have assumed that the SNPs are independent. Below we show
how to specify LD-independent SNPs in the software:

``` r
# Boolean to specify which SNPs are independent (e.g., pruning)
# All SNPs are LD-independent in this example data set 
indep_snps <- rep(TRUE, length(p))

# Create model 
mpi0 <- pi0_model(z = z,
                  indep_snps = indep_snps,
                  knots = c(0.01, 0.025, 0.05, 0.1))

# Estimation fpi0 using design matrix from mpi0
fpi0 <- fpi0est(p = p,
                z = mpi0$zt,
                indep_snps = indep_snps,
                pi0_model = mpi0$fmod)

# Estimate FDR quantities and functional p-value
sffdr_out <- sffdr(p,
                   fpi0 = fpi0$fpi0,
                   indep_snps = indep_snps,
                   epsilon = min(p))
```

Note that the LD-independent SNPs are used for model fitting in sfFDR,
and the functional p-values, q-values, and local FDRs are estimated for
all SNPs. See `?sffdr` for additional details and input arguments and
`?ffinemap` to perform functional fine mapping in a region of interest
(assuming a single causal locus) with the functional local FDR.
