
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Build
Status](https://travis-ci.org/rauschenberger/globalSeq.svg?branch=master)](https://travis-ci.org/rauschenberger/globalSeq)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/rauschenberger/globalSeq?svg=true)](https://ci.appveyor.com/project/rauschenberger/globalSeq)
[![Coverage
Status](https://codecov.io/github/rauschenberger/globalSeq/coverage.svg?branch=master)](https://codecov.io/github/rauschenberger/globalSeq?branch=master)

## Scope

Testing for association between RNA-Seq and other genomic data is
challenging due to high variability of the former and high
dimensionality of the latter.

Using the negative binomial distribution and a random effects model, we
developed an omnibus test that overcomes both difficulties. It may be
conceptualised as a test of overall significance in regression analysis,
where the response variable is overdispersed and the number of
explanatory variables exceeds the sample size.

The proposed test can detect genetic and epigenetic alterations that
affect gene expression. It can examine complex regulatory mechanisms of
gene expression.

## Installation

The package globalSeq depends on [R
\>= 3.0.0](https://cran.r-project.org/), and is available from
[Bioconductor](http://bioconductor.org/packages/globalSeq/):

``` r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("globalSeq")
```

Alternatively, it can be installed from
[GitHub](https://github.com/rauschenberger/globalSeq). This requires the
package
[devtools](https://CRAN.R-project.org/package=devtools):

``` r
devtools::install_github("rauschenberger/globalSeq",build_vignettes=TRUE)
```

Please restart R before loading the package and its documentation:

``` r
library(globalSeq)
utils::help(globalSeq)
utils::vignette("globalSeq")
```

## Reference

A Rauschenberger, MA Jonker, MA van de Wiel, and RX Menezes (2016).
Testing for association between RNA-Seq and high-dimensional data. BMC
Bioinformatics. 17:118.
[html](http://dx.doi.org/10.1186/s12859-016-0961-5)
[pdf](http://www.biomedcentral.com/content/pdf/s12859-016-0961-5.pdf)

[![Platforms](http://www.bioconductor.org/shields/availability/devel/globalSeq.svg)](http://bioconductor.org/packages/devel/bioc/html/globalSeq.html#archives)
[![Downloads](http://www.bioconductor.org/shields/downloads/globalSeq.svg)](http://bioconductor.org/packages/stats/bioc/globalSeq/)
[![Posts](http://www.bioconductor.org/shields/posts/globalSeq.svg)](https://support.bioconductor.org/t/globalseq/)
[![in
Bioc](http://www.bioconductor.org/shields/years-in-bioc/globalSeq.svg)](http://bioconductor.org/packages/devel/bioc/html/globalSeq.html#since)
[![Build](http://www.bioconductor.org/shields/build/devel/bioc/globalSeq.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/globalSeq/)
