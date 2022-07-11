---
author_profile: true
title: "Simulation Study -- Coverage Probability of 95% Confidence Intervals of the Horvitz-Thompson Estimator"
usemathjax: true
toc: true
toc_label: ""
toc_icon: "fas fa-folder-open"
toc_sticky: false
categories:
  - Notes
tags:
  - Survey Sampling
  - Coverage Probability
  - Simulation
---

This short <a href="https://brilliant.org/wiki/monte-carlo/">simulation</a> study 
examines the <a href="https://en.wikipedia.org/wiki/Coverage_probability">coverage probability</a> of 95% confidence intervals 
of the <a href="https://en.wikipedia.org/wiki/Horvitz%E2%80%93Thompson_estimator">Horvitz-Thompson estimator</a> of the population mean 
under simple random sampling without replacement. 

## 95% Confidence Intervals
The Horvitz-Thompson estimator of the population mean:

$$
\begin{align*}
\bar{y} = \frac{1}{n}\sum_{k \in S} y_k
\end{align*}
$$

where $$S$$ is the set of samples and $$n$$ is the sample size.

A corresponding 95% confidence interval:

$$
\begin{align*}
\left( 
  \bar{y} \pm 1.96\sqrt{\hat{Var}[\bar{y}]}
\right)
\end{align*}
$$

where the variance estimator is

$$
\begin{align*}
\hat{Var}[\bar{y}] = \left(1 - \frac{n}{N} \right) \frac{s_y^2}{n}
\end{align*}
$$

and $$s_y^2$$ is the sample variance.

## Simulation
I start with generating the following populations of survey variable $$y$$. Each population is of size N = 10000.

1. N(10, 25)

2. Weibull(shape = 0.5, scale = 4)

The simulation contains R = 50000 replication.

<a href="simulation_function.R">*simulation_function.R*</a> contains a function to execute this simulation.
```R
source("simulation_function.R")
```

Set the population size.
```R
N = 10000
```

### Simulation: N(10, 25)
```R
set.seed(123)
y.population.1 = rnorm(N, mean = 10, sd = 5)
output.1 = simulation(y.population.1, R = 50000)
output.1
```
<table class="table table-bordered table-hover table-condensed">
<thead><tr><th title="Field #1"></th>
<th title="Field #2">sample.size</th>
<th title="Field #3">Rel.Bias.Point.Est</th>
<th title="Field #4">Var.Point.Est</th>
<th title="Field #5">Rel.Bias.Var.Est</th>
<th title="Field #6">Coverage.Prob</th>
</tr></thead>
<tbody><tr>
<td align="right">1</td>
<td>n=10</td>
<td align="right">0.556</td>
<td align="right">25.574</td>
<td align="right">0.27</td>
<td align="right">0.922</td>
</tr>
<tr>
<td align="right">2</td>
<td>n=25</td>
<td align="right">0.116</td>
<td align="right">10.176</td>
<td align="right">-0.104</td>
<td align="right">0.941</td>
</tr>
<tr>
<td align="right">3</td>
<td>n=50</td>
<td align="right">0.104</td>
<td align="right">5.067</td>
<td align="right">-0.264</td>
<td align="right">0.947</td>
</tr>
<tr>
<td align="right">4</td>
<td>n=100</td>
<td align="right">-0.016</td>
<td align="right">2.524</td>
<td align="right">-0.121</td>
<td align="right">0.949</td>
</tr>
<tr>
<td align="right">5</td>
<td>n=500</td>
<td align="right">0.031</td>
<td align="right">0.485</td>
<td align="right">0.042</td>
<td align="right">0.951</td>
</tr>
</tbody></table>

### Simulation: Weibull(shape = 1.5, scale = 1)
```R
set.seed(123)
y.population.2 = rweibull(N, shape = 0.5, scale = 4)
output.2 = simulation(y.population.2, R = 50000)
output.2
```
<table class="table table-bordered table-hover table-condensed">
<thead><tr><th title="Field #1"></th>
<th title="Field #2">sample.size</th>
<th title="Field #3">Rel.Bias.Point.Est</th>
<th title="Field #4">Var.Point.Est</th>
<th title="Field #5">Rel.Bias.Var.Est</th>
<th title="Field #6">Coverage.Prob</th>
</tr></thead>
<tbody><tr>
<td align="right">1</td>
<td>n=10</td>
<td align="right">0.704</td>
<td align="right">31.57</td>
<td align="right">1.839</td>
<td align="right">0.733</td>
</tr>
<tr>
<td align="right">2</td>
<td>n=25</td>
<td align="right">-0.693</td>
<td align="right">12.117</td>
<td align="right">-2.132</td>
<td align="right">0.818</td>
</tr>
<tr>
<td align="right">3</td>
<td>n=50</td>
<td align="right">-0.104</td>
<td align="right">6.151</td>
<td align="right">-0.395</td>
<td align="right">0.867</td>
</tr>
<tr>
<td align="right">4</td>
<td>n=100</td>
<td align="right">-0.095</td>
<td align="right">3.060</td>
<td align="right">-0.402</td>
<td align="right">0.901</td>
</tr>
<tr>
<td align="right">5</td>
<td>n=500</td>
<td align="right">0.095</td>
<td align="right">0.591</td>
<td align="right">0.305</td>
<td align="right">0.938</td>
</tr>
</tbody></table>

## Observations
* As sample size $$n$$ increases, in general, the magnitudes of the relative biases 
  (the second and the fourth columns) of both point estimators $\bar{y}$ and $\hat{Var}[\bar{y}]$ 
  decrease and approach 0. This phenomenon is consistent with law of large numbers. 
  The magnitudes of biases of both estimators shrinks towards zero as we collect more samples.
  
* The point estimates from $$\hat{Var}[\bar{y}]$$ decrease as sample size $$n$$ increases, reflecting 
  the reduction in the typical error of $$\bar{y}$$ as we collect more samples. 
  In other words, $$\bar{y}$$ is more accurate with more samples (note that $\bar{y}$ is unbiased, so variance is equal to MSE). 
  This observation is consistent with the above point.
  
* The coverage probability of the two-sided confidence interval of $$\bar{Y}$$ increases and approaches 0.95 as $$n$$ increases.
  That is, as we collect more samples, in the simulation, the proportion of the confidence intervals computed containing 
  the population mean increases with the proportion approaching 95%. 
  This observation is expected because of a central limit theorem and 1.96 being the 97.5% quantile of the standard normal distribution.
  
* Note that as $$n$$ increases, the width of the confidence interval decreases because of decreasing $$\hat{Var}[\bar{y}]$$. 
  Despite this, the coverage probability increases and approaches 95%. 
  This phenomenon is consistent with the decreasing relative bias of $$\bar{y}$$ and $$\bar{y}$$ getting closer to the population mean.

* When $$n = 10$$, the coverage probability corresponding to the normal population is much higher 
  than that to the Weibull population (0.92 vs. 0.81), because the Weibull population is heavily-skewed 
  while the normal population is much closer to the standard normal.
  
* As $$n$$ increases, a central limit theorem starts to apply and hence coverage probabilities of 
  both population approaches 0.95. Note that when $$n = 50$$, the coverage probability corresponding to the 
  normal population is very close to 0.95 (0.947) while that to the Weibull population is still 0.867. 
  When $$n=500$$, the coverage probability corresponding to the normal population reaches 0.95 while that to 
  the Weibull is still 0.94. Therefore, it takes much more samples for the Weibull distribution to 
  achieve 0.95 through a central limit theorem. 
