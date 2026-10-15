# Methods and Verification Notes

## Current method definition

This project implements a first-order, case-weight Infinitesimal Jackknife (IJ) covariance estimate for the coefficients of an OLS regression. For this calculation, the IJ covariance estimate is algebraically equal to the HC0 sandwich covariance estimator. The script therefore labels the output `IJ_SE` and documents the HC0 identity rather than treating “IJ” and “HC0” as universally interchangeable labels.

## Derivation for the OLS model

Let the OLS estimate be \(\hat\beta=(X^\top X)^{-1}X^\top y\), with residuals \(e_i=y_i-x_i^\top\hat\beta\). Under an infinitesimal perturbation to the case weight for observation \(i\), the first-order derivative of the coefficient estimate at unit weights is

\[
\hat g_i = (X^\top X)^{-1}x_i e_i.
\]

The case-weight IJ covariance estimate used here is

\[
\widehat{\mathrm{Var}}_{IJ}(\hat\beta)
=\sum_{i=1}^{n}\hat g_i\hat g_i^\top.
\]

Substituting the derivative expression gives

\[
\widehat{\mathrm{Var}}_{IJ}(\hat\beta)
=(X^\top X)^{-1}X^\top\operatorname{diag}(e_1^2,\ldots,e_n^2)X(X^\top X)^{-1},
\]

which is the HC0 sandwich covariance form. The script obtains standard errors by taking the square roots of the diagonal elements.

This derivation identifies the exact case-weight IJ convention used in this project. It does not assert that every procedure called an infinitesimal jackknife, or every finite-sample correction, must give this same expression.

## Bootstrap method

The script resamples the 100 observed rows with replacement 1,000 times, refits the same OLS model, and computes the standard deviation of the bootstrap estimates for each coefficient. A fixed random seed supports reproducibility.

## Run verification

The author ran the script in RStudio. The supplied Console output reported:

- Random seed: 354
- Sample size: 100
- Requested bootstrap resamples: 1,000
- Complete bootstrap fits: 1,000 of 1,000
- `results/comparison_table.csv`: `[OK]`
- `results/model_summary.txt`: `[OK]`
- `figures/bootstrap_x1_distribution.png`: `[OK]`

The output table from that run is reproduced in `README.md` and `results/comparison_table.csv`. In this package the output column is named `IJ_SE` to match the definition above. The covariance calculation is unchanged; the naming and documentation make the estimator being used explicit. Run the included script to regenerate the output files.

## Interpretation and limits

- The predictors and response are simulated for demonstration; the analysis is not based on actual student records.
- The true coefficients are known because they were used to generate the response. Estimated coefficients will not generally equal them exactly in a finite sample.
- Standard errors quantify estimated uncertainty; they do not establish causation or substantive importance.
- Conventional OLS, case-weight IJ/HC0, and bootstrap standard errors can differ. The observed differences should be described, not treated as proof that one method is always best.
- Efron (1982), Chapter 6, is cited for the general IJ/influence-function connection. The algebraic OLS identity above is derived explicitly for the case-weight estimator implemented here.

## Reference

Efron, B. (1982). *The Jackknife, the Bootstrap and Other Resampling Plans*. SIAM, Chapter 6, pp. 37–47. https://doi.org/10.1137/1.9781611970319.ch6
