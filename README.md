# The Infinitesimal Jackknife: Variance Estimation in Regression Analysis Using R

**Author:** Zainab Teniola Lukman  
**Course:** STA 354 — Statistical Computing  
**Institution:** University of Ilorin  
**Project type:** Undergraduate statistical computing project

## Project overview

This project uses simulated data to study uncertainty in the coefficients of a multiple linear regression model. The R script fits an ordinary least squares (OLS) model and compares three standard-error estimates:

1. Conventional OLS standard errors.
2. A first-order, case-weight Infinitesimal Jackknife (IJ) estimate.
3. Nonparametric bootstrap standard errors based on 1,000 resamples.

For the OLS model implemented here, the case-weight IJ covariance estimate is calculated by taking the outer products of the observation-level coefficient-derivative contributions. Algebraically, this gives the HC0 sandwich covariance estimator. This is a specific, stated definition for this project; it should not be generalized to every method or finite-sample convention described as an infinitesimal jackknife.

## Research question

How do conventional OLS standard errors, a case-weight IJ estimate, and bootstrap standard errors compare in a simulated multiple linear regression example?

## Simulated data and model

The script generates 100 artificial observations. The predictors are illustrative variables named study hours (`X1`), past exams taken (`X2`), and class attendance (`X3`). They are simulated variables, not real student records.

The data-generating model is:

\[
Y = 2.0 + 3.0X_1 + 1.5X_2 + 2.5X_3 + \varepsilon,
\]

where \(\varepsilon\) is generated from a normal distribution with mean 0 and standard deviation 1. The random seed is 354 so the simulation can be reproduced.

## Methods

### 1. Ordinary least squares

The model `Y ~ X1 + X2 + X3` is fitted with R's `lm()` function. Conventional model-based standard errors are extracted from `summary(model)`.

### 2. Case-weight Infinitesimal Jackknife for OLS

Let \(X\) be the design matrix, \(e_i\) the OLS residual for observation \(i\), and \(x_i\) its design-matrix row written as a column vector. The first-order derivative contribution for case \(i\) is

\[
\hat g_i = (X^\top X)^{-1}x_i e_i.
\]

The script estimates the covariance matrix by summing the outer products of these contributions:

\[
\widehat{\mathrm{Var}}_{IJ}(\hat\beta)
= \sum_{i=1}^{n}\hat g_i\hat g_i^\top.
\]

For this OLS calculation, that expression is algebraically equal to the HC0 sandwich covariance estimator,

\[
(X^\top X)^{-1}X^\top\operatorname{diag}(e_i^2)X(X^\top X)^{-1}.
\]

The reported `IJ_SE` values are the square roots of the covariance matrix's diagonal entries. The relationship is based on the case-weight derivative definition used in this script. Other IJ formulations or finite-sample adjustments may differ.

### 3. Nonparametric bootstrap

The script draws 1,000 samples of size 100 with replacement, refits the model for each sample, and calculates the standard deviation of the bootstrap coefficient estimates. The successful RStudio run reported 1,000 complete bootstrap fits.

## Results from the successful RStudio run

The table below reproduces the comparison printed in the successful RStudio run. The final script labels the influence-based covariance calculation `IJ_SE` and defines it explicitly as the first-order case-weight IJ estimator for this OLS model. The numerical calculation is unchanged from the run that produced these values.

| Coefficient | Estimate | True value | OLS SE | IJ SE | Bootstrap SE |
|---|---:|---:|---:|---:|---:|
| Intercept | 1.19006 | 2.0 | 0.79059 | 0.70687 | 0.72779 |
| X1 | 3.11938 | 3.0 | 0.10120 | 0.08618 | 0.08990 |
| X2 | 1.49744 | 1.5 | 0.04755 | 0.04236 | 0.04484 |
| X3 | 2.51095 | 2.5 | 0.06637 | 0.06263 | 0.06368 |

The estimates are close to the coefficients used to generate the simulated response, though exact equality is not expected in a finite sample. The three standard-error estimates differ because they are calculated using different approaches and assumptions. The comparison alone does not establish that one estimate is always preferable.

## Validation status

The analysis was successfully run in RStudio: all 1,000 bootstrap fits completed, and the output checks confirmed that the comparison table, model summary, and bootstrap histogram were created. The values in the table above are from that run. In this final package, the covariance column is named `IJ_SE` to match the explicitly documented case-weight IJ definition; the numerical calculation itself is unchanged. Run the included script to regenerate the output files.

## Files in this project

- `IJ_analysis.R` — complete analysis script; uses base R and requires no additional packages.
- `METHODS_AND_VERIFICATION.md` — equations, interpretation, and verification notes.
- `RUN_GUIDE.md` — step-by-step instructions for running the script in RStudio.
- `PORTFOLIO_SUMMARY.md` — short portfolio and CV description.
- `MEDIUM_ARTICLE_DRAFT.md` — a readable project write-up.
- `results/` — the script writes `comparison_table.csv` and `model_summary.txt` here.
- `figures/` — the script writes `bootstrap_x1_distribution.png` here.

## Reproduce the analysis

1. Open `IJ_analysis.R` in RStudio.
2. Choose **Session → Set Working Directory → To Source File Location**.
3. Click **Source** to run the full script.
4. Check the Console for `[OK]` beside each expected output file.
5. Open the `results/` and `figures/` folders to inspect the generated files.

## Reference

Efron, B. (1982). *The Jackknife, the Bootstrap and Other Resampling Plans*. Society for Industrial and Applied Mathematics, Chapter 6, “The Infinitesimal Jackknife, the Delta Method and the Influence Function,” pp. 37–47. https://doi.org/10.1137/1.9781611970319.ch6

Efron's chapter describes the general connection between the infinitesimal jackknife and influence-function methods. The OLS covariance identity shown above is derived for the case-weight calculation implemented in this project.
