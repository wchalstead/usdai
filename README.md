
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

The usdai gives methods for dimension agnostic hypothesis testing.

``` r
library(usdai)
```

Consider a situation in which we have $`100`$ realizations of a
multivariate normal vector in $`\mathbb{R}^{10}`$, and we want to test
whether $`\mu = 0`$. This is equivalent to testing that the $`\ell_4`$
norm of the vector $`\mu`$ is zero.

We can perform such a test as follows:

First, we generate data with true mean $`\mu = 0`$

``` r
set.seed(1232)
# Generate Data
data <- matrix(rnorm(200 * 10), 200, 10)
```

And we can use the `crossWStatL4` function to get a test statistic, and
a p-value can be recovered using the `pcrossW` which approximates the
CDF of the asymptotic distribution.

``` r
W.test <- crossWStatL4(data, 100)
pval <- 1 - pcrossW(W.test)
paste(pval)
#> [1] "0.555922"
```

With a large p-value as above, we fail to reject the null hypothesis
that $`\mu = 0`$.

Now consider a situation with true mean $`\mu = 0.2`$

``` r
set.seed(1213)
# Generate Data
data <- matrix(rnorm(200 * 10, mean = 0.2), 200, 10)
```

And we can use the `crossWStatL4` function to get a test statistic, and
a p-value can be recovered using the `pcrossW` which approximates the
CDF of the asymptotic distribution.

``` r
W.test <- crossWStatL4(data, 100)
pval <- 1 - pcrossW(W.test)
paste(pval)
#> [1] "0.015059"
```

With the small p-value as above, we reject the null hypothesis
$`\mu = 0`$.

Beyond this mean testing example above, the usdai package also includes
functions such as `crossWStatL4` that allow for variance testing as
well. Additionally, the functions `qcrossW`, `pcrossW`, `rcrossW`, and
`dcrossW` allow for quantile, distribution and random generation of the
asymptotic testing distribution.

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
