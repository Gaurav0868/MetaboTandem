
# MetaboTandem

<!-- badges: start -->
<!-- badges: end -->

The goal of MetaboTandem is to analyze and annotate LC-MS/MS data

## Installation

You can install the development version of MetaboTandem from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Coayala/MetaboTandem")
```

Conversely you can run MetaboTandem using a docker container with:

``` bash
docker pull coayala/metabotandem
```

## How to start the app

For both the container or locally just

``` r
library(MetaboTandem)

MetaboTandemApp()
```

