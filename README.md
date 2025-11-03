
<!-- README.md is generated from README.Rmd. Please edit that file -->

# usdai

<!-- badges: start -->

<!-- badges: end -->

Typical asymptotic theory for hypothesis testing requires that our
dimension $`d`$ remain constant as the sample size $`n`$ grows. `usdai`
provides methods for conducting hypothesis tests about the mean and
variance of high-dimension datasets in absence of this assumption,
requiring no prior knowledge of the behavior of $`d`$ as $`n`$ goes to
infinity.

For more information on cross U-statistics, see Kim and Ramdas (2024).
For more information on the self-normalization process, see Wang and
Shao (2020) and Shao (2015).

## Installation

You can install the development version of usdai from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("wchalstead/usdai")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(usdai)
## basic example code
```

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-kim2024dimension" class="csl-entry">

Kim, Ilmun, and Aaditya Ramdas. 2024. “Dimension-Agnostic Inference
Using Cross u-Statistics.” *Bernoulli* 30 (1): 683–711.

</div>

<div id="ref-shao2015self" class="csl-entry">

Shao, Xiaofeng. 2015. “Self-Normalization for Time Series: A Review of
Recent Developments.” *Journal of the American Statistical Association*
110 (512): 1797–817.

</div>

<div id="ref-wang2020hypothesis" class="csl-entry">

Wang, Runmin, and Xiaofeng Shao. 2020. “Hypothesis Testing for
High-Dimensional Time Series via Self-Normalization.”

</div>

</div>
