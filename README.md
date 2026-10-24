# DTS 204 Lab Assessment: Palmer Penguins Analysis

Statistical analysis of the Palmer Penguins dataset in R (DTS 204: Statistical
Computing, Inference and Modelling, MIVA Open University).

## What is in this repository

| File / folder | Description |
|:--|:--|
| `penguins_analysis.R` | Complete R script for Tasks 1-5 |
| `DTS204_Lab_Report.pdf` | Full written report (Task 6) |
| `figures/` | Charts produced by the script (PNG) |
| `tables/` | Result tables produced by the script (CSV) |

## What the analysis covers

1. **Data exploration:** structure, missing values, summary statistics
2. **Hypothesis testing:** Welch t-test (Adelie vs Gentoo body mass) and a one-sample proportion test (share of males)
3. **Correlation:** correlation matrix, tests and heatmap
4. **Regression:** simple and multiple linear regression, model comparison and diagnostics
5. **Visualization:** scatter plot with regression line, boxplot, histogram

## How to run

1. Install [R](https://cran.r-project.org) (4.x or later).
2. Open this folder in VS Code or RStudio so it is the working directory.
3. In the R console run:

   ```r
   install.packages(c("palmerpenguins", "ggplot2", "dplyr"))
   source("penguins_analysis.R")
   ```

Results print in the console; figures and tables are saved in `figures/` and `tables/`.
The script uses `palmerpenguins::penguins` explicitly, because recent R versions
include a different built-in dataset also called `penguins`.


