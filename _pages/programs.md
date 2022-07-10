---
author_profile: true
layout: single
title: Some Programs
use_math: true
permalink: /programs/
toc: true
toc_label: "Some Programs"
toc_icon: "fas fa-folder-open"
toc_sticky: false
---

## Sampling Method


## Forecasting
1. <a href="https://github.com/johntwk/Diebold-Mariano-Test">The Diebold-Mariano Test for Forecast Equivalance</a> [Python]
: This Python function `dm_test` implements the Diebold-Mariano Test (1995) with modification suggested by Harvey et. al (1997) to statitsitcally identify forecast accuracy equivalance for 2 sets of predictions.
: Jorge Sandoval, a technical consultant at a Brazilian state government, 
  wrote a <a href="https://www.kaggle.com/code/jorgesandoval/xgboost-vs-lightgbm-using-diebold-mariano-test/notebook">notebook</a> 
  on Kaggle using this code to compare XGBoost and LightGBM.

2. <a href="https://github.com/johntwk/STATA-Backtesting">Backtesting</a> [Stata]
:  A function for backtesting forecasts, accompanied by a STATA help file explaining how to use the program.

3. <a href="https://github.com/johntwk/Python-ML-rolling-grid-search">Grid Search Hyperparameter Tuning and Rolling Forecast</a> [Python]
: This Python function uses machine learning modelling object from scikit-learn to implement a design of grid search hyperparameter selection that respects temporal ordering of time-series, and forecast time-series using the sliding (rolling)-window strategy.

## Model
1. <a href="https://www.dropbox.com/s/jf9vq9rc12hh3xk/TVAR.zip?dl=0">Threshold Vector Autoregression with 2 Regimes</a> [Stata]
: `TVAR_2r_grid_search` implements a grid search of the optimal threshold variable according to a given criterion.
: `TVAR_2r estimates` a TVAR(2) (with the first regime indexed as 0 and the second regime indexed as 1) using the Ordinary Least Square (OLS) Estimation or the Maximum Likelihood (ML) Estimation. 
: A collaborative effort with Myeongwan Kim.

## Others
1. <a href="https://github.com/johntwk/summary-statistics-sparkline">Summary Statistics and Sparkline Generation</a> [Python]
: This Python class takes pandas dataframe as input to generate formatted summary statistics outputs in text, Latex and PDF with sparklines.

2. <a href="https://github.com/johntwk/Capital-Asset-Pricing-Model-CAPM-">Capital Asset Pricing Model</a> [Excel VBA]

3. <a href="https://github.com/johntwk/Improved-String-Class">An Improved String Class</a> [C++]
: This C++ class allows the overloading of multiplication, accompanied by a set of test cases.


