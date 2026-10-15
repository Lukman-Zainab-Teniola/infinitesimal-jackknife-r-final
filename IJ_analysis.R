# ============================================================
# The Infinitesimal Jackknife: Variance Estimation in Regression
# Undergraduate Statistical Computing project | STA 354
# University of Ilorin
#
# METHOD NOTE:
# This script estimates coefficient uncertainty using the first-order
# case-weight Infinitesimal Jackknife (IJ) for OLS. For this OLS model,
# summing the outer products of the case-weight derivative contributions
# gives the HC0 sandwich covariance matrix. This equivalence is specific
# to the calculation defined here; other IJ procedures or finite-sample
# corrections may use different conventions.
# ============================================================

# 1. Setup -----------------------------------------------------
rm(list = ls())
seed <- 354
set.seed(seed)

n <- 100
B <- 1000

# In RStudio, first choose:
# Session > Set Working Directory > To Source File Location
# This makes results/ and figures/ appear beside this script.
cat("Working directory:", getwd(), "\n")

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

# 2. Simulate data --------------------------------------------
X1 <- rnorm(n, mean = 5, sd = 1)       # Study hours
X2 <- rnorm(n, mean = 3, sd = 2)       # Past exams taken
X3 <- rnorm(n, mean = 8, sd = 1.5)     # Class attendance

beta0 <- 2.0
beta1 <- 3.0
beta2 <- 1.5
beta3 <- 2.5

epsilon <- rnorm(n, mean = 0, sd = 1)
Y <- beta0 + beta1 * X1 + beta2 * X2 + beta3 * X3 + epsilon
dat <- data.frame(Y = Y, X1 = X1, X2 = X2, X3 = X3)

# 3. Ordinary least squares regression ------------------------
model <- lm(Y ~ X1 + X2 + X3, data = dat)
model_summary <- summary(model)
beta_hat <- coef(model)
residuals_ols <- residuals(model)
ols_se <- coef(model_summary)[, "Std. Error"]

# 4. Case-weight Infinitesimal Jackknife covariance -------------
X_matrix <- model.matrix(model)
XtX_inv <- solve(crossprod(X_matrix))

# For observation i, the derivative of the OLS coefficient vector with
# respect to its case weight at unit weights is:
#   g_i = (X'X)^(-1) x_i e_i
# The rows below contain g_i' (equivalently e_i x_i' (X'X)^(-1)).
# The case-weight IJ covariance estimate is sum_i g_i g_i', which is
# algebraically the HC0 sandwich covariance estimator for this OLS model.
case_weight_derivatives <- X_matrix * as.numeric(residuals_ols)
case_weight_derivatives <- case_weight_derivatives %*% XtX_inv
ij_covariance <- crossprod(case_weight_derivatives)
ij_se <- sqrt(pmax(diag(ij_covariance), 0))

# 5. Bootstrap standard errors --------------------------------
coefficient_names <- names(coef(model))
boot_coefs <- matrix(
  NA_real_, nrow = B, ncol = length(coefficient_names),
  dimnames = list(NULL, coefficient_names)
)

for (b in seq_len(B)) {
  indices <- sample(seq_len(n), size = n, replace = TRUE)
  boot_data <- dat[indices, , drop = FALSE]
  boot_model <- try(lm(Y ~ X1 + X2 + X3, data = boot_data), silent = TRUE)

  if (!inherits(boot_model, "try-error")) {
    boot_estimates <- coef(boot_model)
    boot_coefs[b, names(boot_estimates)] <- boot_estimates
  }
}

valid_bootstrap_rows <- complete.cases(boot_coefs)
number_valid_bootstraps <- sum(valid_bootstrap_rows)
if (number_valid_bootstraps < 2) {
  stop("Fewer than two complete bootstrap fits; bootstrap SEs cannot be estimated.")
}
boot_se <- apply(boot_coefs, 2, sd, na.rm = TRUE)

# 6. Comparison table and text report -------------------------
comparison <- data.frame(
  Coefficient = coefficient_names,
  Estimate = as.numeric(beta_hat[coefficient_names]),
  True_Value = c(beta0, beta1, beta2, beta3),
  OLS_SE = as.numeric(ols_se[coefficient_names]),
  IJ_SE = as.numeric(ij_se[coefficient_names]),
  Bootstrap_SE = as.numeric(boot_se[coefficient_names]),
  row.names = NULL
)
numeric_columns <- vapply(comparison, is.numeric, logical(1))
comparison[numeric_columns] <- lapply(comparison[numeric_columns], round, 5)

write.csv(comparison, file.path("results", "comparison_table.csv"), row.names = FALSE)

capture.output({
  cat("REPRODUCIBLE REGRESSION ANALYSIS\n\n")
  cat("Random seed:", seed, "\n")
  cat("Sample size:", n, "\n")
  cat("Requested bootstrap resamples:", B, "\n")
  cat("Complete bootstrap fits:", number_valid_bootstraps, "\n\n")
  cat("OLS MODEL SUMMARY\n")
  print(model_summary)
  cat("\nCOMPARISON TABLE\n")
  print(comparison, row.names = FALSE)
  cat("\nMethod note:\n")
  cat(
    "IJ_SE is based on the first-order case-weight derivative contributions ",
    "for OLS. Summing their outer products gives the HC0 sandwich ",
    "covariance estimate for this model. This statement refers to the ",
    "case-weight IJ definition implemented in this script; other IJ ",
    "procedures or finite-sample conventions may differ.\n", sep = ""
  )
}, file = file.path("results", "model_summary.txt"))

# 7. Bootstrap histogram for X1 -------------------------------
plot_file <- file.path("figures", "bootstrap_x1_distribution.png")
x1_boot <- boot_coefs[, "X1"]
x1_boot <- x1_boot[is.finite(x1_boot)]
if (length(x1_boot) < 2) {
  stop("Not enough valid X1 bootstrap estimates to draw the histogram.")
}

tryCatch({
  png(filename = plot_file, width = 1200, height = 800, res = 140)
  hist(
    x1_boot,
    main = "Bootstrap Distribution of the Study-Hours Coefficient",
    xlab = "Bootstrap coefficient estimate for X1",
    col = "lightblue",
    border = "white"
  )
  abline(v = beta1, lwd = 2, lty = 2, col = "red")
  abline(v = beta_hat["X1"], lwd = 2, col = "blue")
  legend(
    "topright",
    legend = c("True coefficient", "OLS estimate"),
    col = c("red", "blue"), lwd = 2, lty = c(2, 1), bty = "n"
  )
  dev.off()
}, error = function(e) {
  if (grDevices::dev.cur() > 1) try(grDevices::dev.off(), silent = TRUE)
  stop(
    "Could not create ", plot_file,
    ". Check the working directory and folder permissions. Error: ",
    conditionMessage(e)
  )
})

# 8. Confirm files exist --------------------------------------
expected_files <- c(
  file.path("results", "comparison_table.csv"),
  file.path("results", "model_summary.txt"),
  plot_file
)
cat("\nAnalysis finished.\n")
cat("Working directory:", normalizePath(getwd()), "\n")
cat("Output file checks:\n")
for (f in expected_files) {
  cat(if (file.exists(f)) "[OK] " else "[MISSING] ", f, "\n")
}
cat("\nComplete bootstrap fits:", number_valid_bootstraps, "of", B, "\n\n")
print(comparison, row.names = FALSE)
