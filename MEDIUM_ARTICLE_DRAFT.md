# Comparing Regression Standard Errors in R: An Undergraduate Simulation Project

*By Zainab Teniola Lukman*

## Why uncertainty matters

A regression coefficient is estimated from data, so it is useful to measure how much that estimate may vary. Standard errors are one way to summarize estimated uncertainty in regression coefficients.

For my STA 354 Statistical Computing project at the University of Ilorin, I wrote an R script to compare three approaches in a simulated multiple-regression example: conventional OLS standard errors, a case-weight Infinitesimal Jackknife estimate, and bootstrap standard errors.

## The simulation

The script generates 100 artificial observations with three predictors: study hours (`X1`), past exams taken (`X2`), and class attendance (`X3`). The outcome follows a specified linear model plus normally distributed random error. These variables are simulated and do not represent actual student records. The random seed is set to 354 for reproducibility.

## The three estimates

The first method is the conventional model-based standard error from `lm()`. The second is a first-order case-weight Infinitesimal Jackknife (IJ) covariance estimate, obtained by summing the outer products of observation-level derivatives of the regression coefficients. For the OLS model used here, this expression is algebraically equal to the HC0 sandwich covariance estimator. The third method uses 1,000 nonparametric bootstrap samples and calculates the standard deviation of the refitted coefficients.

This equivalence is specific to the case-weight IJ calculation implemented in this project. It should not be taken to mean that every procedure called an infinitesimal jackknife uses the same formula.

## Results

The successful RStudio run completed all 1,000 bootstrap fits and generated the comparison table, model summary, and bootstrap histogram. The table in the README records the coefficient estimates and standard errors printed by that run. The estimated coefficients were close to the known data-generating values, although exact equality is not expected in a finite sample. The three standard-error approaches produced somewhat different estimates of uncertainty.

## What the project demonstrates

The project links statistical ideas to a reproducible programming workflow: simulate data, fit a model, estimate coefficient uncertainty, resample observations, export results, and visualize a bootstrap distribution. It also illustrates why statistical method labels should be tied to a clearly stated formula.

## Reference

Efron, B. (1982). *The Jackknife, the Bootstrap and Other Resampling Plans*. SIAM, Chapter 6, “The Infinitesimal Jackknife, the Delta Method and the Influence Function.” https://doi.org/10.1137/1.9781611970319.ch6
