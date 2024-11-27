
<!-- README.md is generated from README.Rmd. Please edit that file -->

# EGHS

<!-- badges: start -->
<!-- badges: end -->

The goal of EGHS is to provide two working Harvest Control Rules to
illustrate the use of the aMSE software.

Harvest Control Rules differ from harvest strategies in that they
literally take the stock status, as derived from whatever analysis is
done to the data from the fishery (in the case of abalone this relates
to standardizing the CPUE) and in their turn, derive a future catch
level. A harvest strategy would include descriptions of what data to
use, what analysis to apply, and what harvest control rule to apply. The
current Tasmanian harvest strategy differs from the two harvest control
rules (hcr) included in **EGHS**, and continues to be developed. Hence
it is not suitable for inclusion in the examples within *aMSEGuide*.

Two harvest control rules are provided:

- Constant Catch *consthcr*, which, as its name suggests simply applies
  a constant aspirational catch to each sau across all projections. This
  is the simplest hcr.
- Adaptive reference period, *mcdahcr*, which reflects the original HCR
  used in Tasmania (Bradshaw, 2018), where the reference period used to
  generate CPUE targets increases as the years of projection increase.
  This is a much more complex empirical harvest strategy, especially
  when including the three meta-rules that are used to tune the outcome
  of the hcr to more closely match the objectives of the harvest
  strategy (HS).

## Installation

You can install the development version of EGHS from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("haddonm/EGHS")
```

## Example

``` r
library(EGHS)
?EGHS
#> No documentation for 'EGHS' in specified packages and libraries:
#> you could try '??EGHS'
```