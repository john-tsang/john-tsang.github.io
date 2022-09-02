---
author_profile: true
published: true
title: "Accelerate `R` with `Rcpp`: Sampling with Replacement and the Hansen-Hurwitz Estimator [Draft]"
use_math: true
toc: true
toc_label: ""
toc_icon: "fas fa-folder-open"
toc_sticky: false
categories:
  - Notes
tags:
  - R
  - Rcpp
  - C++
  - OpenMP
  - Survey Sampling
---

This note aims to demonstrate the use of `Rcpp` and `OpenMP` to accelerate computation in `R`.

# `Rcpp` for Performance in `R`
[`Rcpp`](https://www.rcpp.org/) is an `R` package that allows for using the more efficient yet more 
complicated [`C++`](https://en.wikipedia.org/wiki/C%2B%2B) with `R` to improve computational time.
The improvement is usually tremendous, so it is worthwhile to consider [`Rcpp`](https://www.rcpp.org/) for large computational tasks.

# Task 1: Sampling with Replacement
Consider the case of sampling with replacement (equal probability of being chosen) from a sample of size 100,000. 

Generate sampling frame:
```R
set.seed(123465)
n.elem = 100000
frame1 = rnorm(n.elem)
```
## R Vs. Rcpp: Rcpp Faster
`sample_with_replacement.cpp`: 
```c++
#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export()]]
arma::vec sample_w_replacement(arma::vec x, const int size) {
  //int size = x.n_rows;
  arma::uvec index = arma::randi<arma::uvec>(size, arma::distr_param(0,size-1));
  arma::vec rt = x.elem(index);
  return rt;
}
```

Load the library `Rcpp` and the `C++` source file.
```R
library(Rcpp)
sourceCpp("sample_with_replacement.cpp")
```


Runtime comparison: `sample` in `R` Vs. `sample_w_replacement` with `Rcpp`
```R
library(microbenchmark)
microbenchmark(
  r = sample(frame1, n.elem, replace = TRUE),
  cpp = sample_w_replacement(frame1, n.elem),
  times = 1e3
)
```

*Output:*
```
Unit: milliseconds
 expr    min      lq     mean  median     uq     max neval
    r 7.2487 8.14530 9.007408 8.63335 9.3552 34.6576  1000
  cpp 1.3987 1.62835 2.294706 1.92230 2.5577 23.6820  1000
```

The `Rcpp` implementation of sampling with replacement runs faster than `sample` in `R`.

# Task 2: A Simulation Study for the Hansen-Hurwitz Estimator

## The Hansen-Hurwitz (HH) Estimator

```R
HH.estimator = function(sample.selected) {
  return(sum(sample.selected) / length(sample.selected))
}
```

## `R` Vs. `Rcpp` Vs. `Rcpp` & `OpenMP`
`R`:
```R
r.sim = function(R, n, frame1) {
  results = 0
  for (i in 1:R) {
    sample.selected = sample(frame1, n, replace = TRUE)
    results = results + HH.estimator(sample.selected)
  }
  abs.rel.bias = (results/R - mean(frame1)) / mean(frame1) * 100
  abs.rel.bias
}
```

`C++`:
```c++
double HH_estimator(arma::vec sample) {
  return(arma::mean(sample));
}

// [[Rcpp::export()]]
double rcpp_sim(const int R, const int n, arma::vec frame1) {
  double results = 0;
  for (int i = 0; i < R; ++i) {
    arma::vec sample_selected = sample_w_replacement(frame1, n);
    results += HH_estimator(sample_selected);
  }
  double abs_rel_bias = (results/R - arma::mean(frame1)) / arma::mean(frame1) * 100;
  return abs_rel_bias;
}
```
Runtime comparison:
```R
microbenchmark(
  r = r.sim(1000, 10000, frame1),
  cpp = rcpp_sim(1000, 10000, frame1),
  times = 1000L
)
```

*Output:*
```
Unit: milliseconds
         expr      min        lq      mean    median       uq
            r 794.5302 1001.8973 1077.5736 1051.0002 1148.692
          cpp 183.4560  219.5191  248.2792  240.3830  269.853
       max neval
 1398.4725  1000
  469.8277  1000
```
The `Rcpp` implementation runs faster than `sample` in `R`.