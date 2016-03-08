<!-- README.md is generated from README.Rmd. Please edit that file -->
**Scope**

Testing for association between RNA-Seq and other genomic data is challenging due to high variability of the former and high dimensionality of the latter.

Using the negative binomial distribution and a random effects model, we developed an omnibus test that overcomes both difficulties. It may be conceptualised as a test of overall significance in regression analysis, where the response variable is overdispersed and the number of explanatory variables exceeds the sample size.

The proposed test can detect genetic and epigenetic alterations that affect gene expression. It can examine complex regulatory mechanisms of gene expression.

**Installation**

The package globalSeq depends on [R \>= 3.3](https://cran.r-project.org/), and is available from [Bioconductor](http://bioconductor.org/packages/globalSeq/):

``` r
source("https://bioconductor.org/biocLite.R")
BiocInstaller::biocLite("globalSeq")
```

Alternatively, it can be installed from [GitHub](https://github.com/Bioconductor-mirror/globalSeq). This requires the package [devtools](https://cran.r-project.org/web/packages/devtools/README.html):

``` r
devtools::install_github("Bioconductor-mirror/globalSeq")
```

Please restart R before loading the package and its documentation:

``` r
library(globalSeq)
utils::help(globalSeq)
utils::vignette("globalSeq")
```

**Reference**

A Rauschenberger, MA Jonker, MA van de Wiel, and RX Menezes (2016). Testing for association between RNA-Seq and high-dimensional data. BMC Bioninformatics. 17:118. [html](http://dx.doi.org/10.1186/s12859-016-0961-5) [pdf](http://www.biomedcentral.com/content/pdf/s12859-016-0961-5.pdf)
