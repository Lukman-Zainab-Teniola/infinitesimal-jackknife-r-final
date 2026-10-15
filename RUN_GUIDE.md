# Run Guide: RStudio

This project uses base R and does not require additional packages.

## Run the analysis

1. Open `IJ_analysis.R` in RStudio.
2. Choose **Session → Set Working Directory → To Source File Location**.
3. Click **Source** to run the entire script.
4. Wait for `Analysis finished.` in the Console.
5. Confirm that all three output checks show `[OK]`.
6. Open `results/comparison_table.csv`, `results/model_summary.txt`, and `figures/bootstrap_x1_distribution.png` from the project folder.

## What the script saves

- `results/comparison_table.csv` — coefficient estimates and standard errors.
- `results/model_summary.txt` — model summary and run information.
- `figures/bootstrap_x1_distribution.png` — bootstrap histogram for the `X1` coefficient.

## If an output is missing

Check the full red error message in the Console, run `getwd()` to check the working directory, select **Session → Set Working Directory → To Source File Location** again, and run the entire script using **Source**.

## Expected run information

With the current script settings, the run uses seed 354, 100 observations, and 1,000 bootstrap resamples. Bootstrap estimates can occasionally fail in other datasets; the script reports the number of complete fits.
