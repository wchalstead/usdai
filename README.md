
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

This testing procedure is based on the self-normalized cross
U-statistic— see Wang and Shao (2020) and Shao (2015) for descriptions
of the self-normalization process and Kim and Ramdas (2024) for
description of the cross U-statistic— whose asymptotic density is
plotted in Lobato (2001).

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

## To-Do

I would like to move my current code to C++ to improve the computation
time and add some basic tests to make sure that my computations are
still correct as I transfer my code over. A final thing that I would
like to add is a function to simulate random variables from the
asymptotic distribution, but I will likely need to do some testing to
decide if it is more efficient to try and generate these by generating a
uniform random variable on (0, 1) and passing this as a function to some
pre-computed set of quantiles or finding a way to efficiently integrate
over a random process.

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-kim2024dimension" class="csl-entry">

Kim, Ilmun, and Aaditya Ramdas. 2024. “Dimension-Agnostic Inference
Using Cross u-Statistics.” *Bernoulli* 30 (1): 683–711.

</div>

<div id="ref-lobato2001testing" class="csl-entry">

Lobato, Ignacio N. 2001. “Testing That a Dependent Process Is
Uncorrelated.” *Journal of the American Statistical Association* 96
(455): 1066–76.

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
