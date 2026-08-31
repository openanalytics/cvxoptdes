
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cvxoptdes <img src='man/figures/logo.svg' align="right" alt="cvxoptdes sticker" width="20%"/>

An R-package for *convex optimal experiment design*, using exclusively
open-source optimization solvers. `cvxoptdes` provides an
object-oriented interface (via [R6](https://r6.r-lib.org/)) for
specifying and optimizing designs for linear models, generalized
linear/additive models, nonlinear models, and random intercept linear
mixed models. It supports optimization and augmentation of weighted,
penalized and convex constrained approximate designs across a discrete
grid of candidate points by solving a second-order cone program using
the [Splitting Cone Solver (SCS)](https://www.cvxgrp.org/scs/) and
optimal rounding to integer-valued designs using point-exchange
algorithms or mixed-integer programming with the
[HiGHS](https://highs.dev/) solver for orthogonally constrained designs.

## Installation

When installing the R-package from source, verify that
[CMake](https://cmake.org/) is installed on the system and is found by R
using e.g.,

``` r
Sys.which("cmake")
#>           cmake 
#> "/usr/bin/cmake"
```

With [Cmake](https://cmake.org/) available, install the R-package from
source with:

``` r
## install.packages("remotes")
remotes::install_github("openanalytics/cvxoptdes")
```

### Windows

Using [Rtools42](https://cran.r-project.org/bin/windows/Rtools/) or
higher, [Cmake](https://cmake.org/) comes with the Rtools installation
and is available by default.

## Quick start

Optimal designs for linear models are calculated using the
[R6](https://r6.r-lib.org/articles/Introduction.html)-class `lm_design`.
When initializing a design object with `lm_design$new()`, the minimal
required inputs are a (discrete) candidate design grid and a one-sided
linear model formula:

``` r
library(cvxoptdes)

## 3 factors w/ 5, 5 and 3 levels
candidates <- expand.grid(
  x1 = -2:2,
  x2 = -2:2,
  x3 = -1:1
)

design <- lm_design$new(
  formula = ~(x1 + x2 + x3)^2 + I(x1^2) + I(x2^2) + I(x3^2),   ## full quadratic model
  data = candidates 
)
```

After initializing the linear model and design grid, use the
`$optimize()` method to calculate an *approximate* optimal design:

``` r
design$optimize(criterion = "D")
#>  [1] 0.06 0.00 0.03 0.00 0.06 0.00 0.00 0.00 0.00 0.00 0.03 0.00 0.00 0.00 0.03 0.00 0.00 0.00 0.00 0.00 0.06 0.00 0.03 0.00 0.06 0.03 0.00 0.00 0.00
#> [30] 0.03 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.06 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.03 0.00 0.00 0.00 0.03 0.06 0.00 0.03 0.00 0.06 0.00 0.00 0.00
#> [59] 0.00 0.00 0.03 0.00 0.00 0.00 0.03 0.00 0.00 0.00 0.00 0.00 0.06 0.00 0.03 0.00 0.06
#> attr(,"criterion")
#> attr(,"criterion")$crit.value
#> [1] 1.897913
#> 
#> attr(,"criterion")$crit
#> [1] "D"
```

The `$optimize()` method returns the weights assigned to each design
point in `design$data`, as well as the optimized value of the selected
design criterion.

The `$round()` method then rounds the approximate optimal design to an
exact (integer) design. If `method = "optimal"`, repeated point-exchange
procedures are executed to find the optimal exact design:

``` r
## parallel evaluation
mirai::daemons(n = parallel::detectCores() / 2, seed = 1)   

exact_design <- design$round(m = 13, method = "optimal", n_repeats = 100)

mirai::daemons(0)

exact_design
#>  [1] 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 1 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0
#> [74] 0 1
#> attr(,"criterion")
#> attr(,"criterion")$crit
#> [1] "D"
#> 
#> attr(,"criterion")$crit.value
#> [1] 1.84352
#> 
#> attr(,"criterion")$crit.rel.eff
#> [1] 0.9713409

design$subset(replicates = exact_design)
#>    x1 x2 x3
#> 1  -2 -2 -1
#> 5   2 -2 -1
#> 13  0  0 -1
#> 21 -2  2 -1
#> 25  2  2 -1
#> 28  0 -2  0
#> 36 -2  0  0
#> 50  2  2  0
#> 51 -2 -2  1
#> 55  2 -2  1
#> 65  2  0  1
#> 71 -2  2  1
#> 73  0  2  1
```

The existing `lm_design` object can be reused to compute optimal designs
for other design criteria, or to evaluate other properties of a given
design. For instance, the variance-covariance matrix of the linear model
parameter estimates:

``` r
design$vcov(design_weights = exact_design, sigma2 = 1)
#>             (Intercept)     x1     x2     x3 I(x1^2) I(x2^2) I(x3^2)  x1:x2  x1:x3  x2:x3
#> (Intercept)       0.722 -0.002  0.002  0.005  -0.070  -0.070  -0.279  0.012  0.024 -0.024
#> x1               -0.002  0.025  0.000  0.000   0.001   0.000   0.000  0.000 -0.001  0.003
#> x2                0.002  0.000  0.025  0.000   0.000  -0.001   0.000  0.000  0.003 -0.001
#> x3                0.005  0.000  0.000  0.101   0.000   0.000  -0.006  0.003  0.001 -0.001
#> I(x1^2)          -0.070  0.001  0.000  0.000   0.030  -0.003  -0.014  0.000  0.000  0.008
#> I(x2^2)          -0.070  0.000 -0.001  0.000  -0.003   0.030  -0.014  0.000 -0.008  0.000
#> I(x3^2)          -0.279  0.000  0.000 -0.006  -0.014  -0.014   0.474 -0.015  0.000  0.000
#> x1:x2             0.012  0.000  0.000  0.003   0.000   0.000  -0.015  0.009  0.003 -0.003
#> x1:x3             0.024 -0.001  0.003  0.001   0.000  -0.008   0.000  0.003  0.035 -0.005
#> x2:x3            -0.024  0.003 -0.001 -0.001   0.008   0.000   0.000 -0.003 -0.005  0.035
```

## Links

Learn more at <https://repos.openanalytics.eu/docs/cvxoptdes/>

------------------------------------------------------------------------

**(c) Copyright Open Analytics NV, 2026 - Apache License 2.0**
