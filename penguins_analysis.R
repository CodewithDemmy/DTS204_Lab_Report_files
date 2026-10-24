# =============================================================================
# DTS 204 - Statistical Computing, Inference and Modelling
# LAB ASSESSMENT: Statistical analysis of the Palmer Penguins dataset
#
# Author      : Ademola Utibe Afolabi (Demmy)
# Matric No.  : 2025/A/CSC/0096
# Institution : MIVA Open University
#
# HOW TO RUN
#   1. Open the folder containing this file in VS Code / RStudio.
#   2. Make sure the working directory is that folder
#      (VS Code: open the folder, then run in the R terminal;
#       RStudio: Session > Set Working Directory > To Source File Location).
#   3. Run the whole script:   source("penguins_analysis.R")
#      (or run it section by section, top to bottom).
#
# OUTPUT
#   - Printed results in the console (hypotheses, p-values, models, ...)
#   - figures/  : all plots as PNG files
#   - tables/   : key result tables as CSV files
#
# CONTENTS
#   0. Setup
#   1. Task 1: Data exploration and preprocessing
#   2. Task 2: Hypothesis testing (t-test, proportion test)
#   3. Task 3: Correlation analysis
#   4. Task 4: Regression modelling (simple and multiple)
#   5. Task 5: Visualization
#   6. Session information
# =============================================================================


# =============================================================================
# 0. SETUP
# =============================================================================

# Install any missing packages (needs an internet connection the first time).
required_pkgs <- c("palmerpenguins", "ggplot2", "dplyr")
missing_pkgs  <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(palmerpenguins)
  library(ggplot2)
  library(dplyr)
})

# Folders for outputs
dir.create("figures", showWarnings = FALSE)
dir.create("tables",  showWarnings = FALSE)

# Significance level used for every test in this analysis
alpha <- 0.05

# Colour-blind-friendly colours for the three species (used in all plots)
species_cols <- c(Adelie = "#E69F00", Chinstrap = "#0072B2", Gentoo = "#009E73")

# Small helper to print clear section headings in the console
section <- function(title) {
  cat("\n", strrep("=", 70), "\n", title, "\n", strrep("=", 70), "\n", sep = "")
}

# Helper: format a p-value (scientific notation when very small)
fmt_p <- function(p) {
  if (p < 0.001) formatC(p, format = "e", digits = 2)
  else formatC(p, format = "f", digits = 4)
}


# =============================================================================
# 1. TASK 1: DATA EXPLORATION AND PREPROCESSING
# =============================================================================
section("TASK 1: DATA EXPLORATION AND PREPROCESSING")

# --- 1.1 Load the dataset ----------------------------------------------------
# NOTE: palmerpenguins:: is written out in full on purpose. Recent R versions
# (4.5 and later) also include a *different* dataset called `penguins` in base
# R (with column names such as bill_len), so we explicitly ask for the
# palmerpenguins version, whose columns match the lab brief.
penguins <- palmerpenguins::penguins

cat("\nDimensions (rows, columns):", dim(penguins), "\n\n")
cat("Structure of the data:\n")
str(penguins)

# --- 1.2 First six rows ------------------------------------------------------
cat("\nFirst six rows:\n")
print(head(penguins, 6))

# --- 1.3 Missing values ------------------------------------------------------
missing_summary <- data.frame(
  variable    = names(penguins),
  n_missing   = as.integer(colSums(is.na(penguins))),
  pct_missing = round(100 * colMeans(is.na(penguins)), 2)
)
cat("\nMissing values per variable:\n")
print(missing_summary)

incomplete_rows <- penguins[!complete.cases(penguins), ]
cat("\nNumber of rows with at least one missing value:",
    nrow(incomplete_rows), "of", nrow(penguins), "\n")
cat("Number of fully complete rows:", sum(complete.cases(penguins)), "\n\n")
cat("Rows with missing values:\n")
print(incomplete_rows, n = Inf)

# How the missing values are handled
# ---------------------------------------------------------------------------
# * Two penguins have NO measurements at all (bill length, bill depth,
#   flipper length and body mass are all NA). They carry no information for
#   the numeric analyses and are excluded from them.
# * Nine further penguins have a measured body but an unknown `sex`.
# * The amount of missing data is very small (under 1% for the measurements,
#   about 3% for sex), so I use ANALYSIS-SPECIFIC (available-case) deletion:
#   each test / model uses every row that has the variables it needs.
#   - Measurement analyses (summary stats, correlations, regressions, t-test)
#     use the 342 penguins with measurements.
#   - The proportion test for sex uses the 333 penguins with a recorded sex.
# * I do NOT impute sex: guessing it would bias the proportion test, which is
#   exactly the quantity being tested.

# --- 1.4 Summary statistics for numerical variables -------------------------
num_vars <- c("bill_length_mm", "bill_depth_mm",
              "flipper_length_mm", "body_mass_g")

summary_stats <- data.frame(
  variable = num_vars,
  n        = sapply(penguins[num_vars], function(x) sum(!is.na(x))),
  mean     = sapply(penguins[num_vars], mean,   na.rm = TRUE),
  median   = sapply(penguins[num_vars], median, na.rm = TRUE),
  min      = sapply(penguins[num_vars], min,    na.rm = TRUE),
  max      = sapply(penguins[num_vars], max,    na.rm = TRUE),
  sd       = sapply(penguins[num_vars], sd,     na.rm = TRUE),
  row.names = NULL
)
summary_stats$range <- summary_stats$max - summary_stats$min
summary_stats <- summary_stats[, c("variable", "n", "mean", "median",
                                   "min", "max", "range", "sd")]
summary_stats[-1] <- round(summary_stats[-1], 2)

cat("\nSummary statistics (mean, median, range, standard deviation):\n")
print(summary_stats)

# Extra detail: body mass broken down by species
body_mass_by_species <- penguins %>%
  group_by(species) %>%
  summarise(
    n      = sum(!is.na(body_mass_g)),
    mean   = round(mean(body_mass_g,   na.rm = TRUE), 1),
    median = round(median(body_mass_g, na.rm = TRUE), 1),
    sd     = round(sd(body_mass_g,     na.rm = TRUE), 1),
    min    = min(body_mass_g, na.rm = TRUE),
    max    = max(body_mass_g, na.rm = TRUE),
    .groups = "drop"
  )
cat("\nBody mass (g) by species:\n")
print(as.data.frame(body_mass_by_species))

write.csv(missing_summary,       "tables/t1_missing_values.csv",   row.names = FALSE)
write.csv(summary_stats,         "tables/t1_summary_stats.csv",    row.names = FALSE)
write.csv(body_mass_by_species,  "tables/t1_body_mass_species.csv", row.names = FALSE)


# =============================================================================
# 2. TASK 2: HYPOTHESIS TESTING
# =============================================================================
section("TASK 2: HYPOTHESIS TESTING")

# --- 2.i Two-sample t-test: body mass of Adelie vs Gentoo --------------------
# H0: mean body mass (Adelie) = mean body mass (Gentoo)
# H1: mean body mass (Adelie) != mean body mass (Gentoo)   (two-sided)
cat("\n--- 2.i Two-sample t-test: body mass, Adelie vs Gentoo ---\n")

adelie <- penguins$body_mass_g[penguins$species == "Adelie" &
                                 !is.na(penguins$body_mass_g)]
gentoo <- penguins$body_mass_g[penguins$species == "Gentoo" &
                                 !is.na(penguins$body_mass_g)]

cat("\nSample sizes  : Adelie =", length(adelie), " Gentoo =", length(gentoo))
cat("\nSample means  : Adelie =", round(mean(adelie), 1),
    " Gentoo =", round(mean(gentoo), 1))
cat("\nSample SDs    : Adelie =", round(sd(adelie), 1),
    " Gentoo =", round(sd(gentoo), 1), "\n")

# Assumption checks
cat("\nAssumption check 1 - normality within each group (Shapiro-Wilk):\n")
shapiro_adelie <- shapiro.test(adelie)
shapiro_gentoo <- shapiro.test(gentoo)
print(shapiro_adelie)
print(shapiro_gentoo)

cat("\nAssumption check 2 - equal variances (F test):\n")
var_res <- var.test(adelie, gentoo)
print(var_res)

# Main test: Welch two-sample t-test (does not assume equal variances)
t_res <- t.test(adelie, gentoo, alternative = "two.sided",
                var.equal = FALSE, conf.level = 1 - alpha)
cat("\nMAIN TEST - Welch two-sample t-test:\n")
print(t_res)

mean_diff <- unname(t_res$estimate[1] - t_res$estimate[2])
pooled_sd <- sqrt(((length(adelie) - 1) * var(adelie) +
                     (length(gentoo) - 1) * var(gentoo)) /
                    (length(adelie) + length(gentoo) - 2))
cohens_d  <- mean_diff / pooled_sd

cat("\nMean difference (Adelie - Gentoo):", round(mean_diff, 2), "g\n")
cat("95% CI for the difference        :",
    round(t_res$conf.int[1], 2), "to", round(t_res$conf.int[2], 2), "g\n")
cat("Cohen's d (effect size)          :", round(cohens_d, 2), "\n")
cat("p-value                          :", fmt_p(t_res$p.value), "\n")
cat("Decision at alpha =", alpha, ":",
    ifelse(t_res$p.value < alpha,
           "REJECT H0 - mean body mass differs significantly.",
           "FAIL TO REJECT H0 - no significant difference."), "\n")

# Robustness checks (same conclusion expected)
cat("\nRobustness check A - Student t-test (equal variances assumed):\n")
t_student <- t.test(adelie, gentoo, var.equal = TRUE)
print(t_student)

cat("\nRobustness check B - Wilcoxon rank-sum test (no normality assumption):\n")
wilcox_res <- suppressWarnings(wilcox.test(adelie, gentoo))
print(wilcox_res)

t_test_table <- data.frame(
  test      = c("Welch t-test", "Student t-test", "Wilcoxon rank-sum"),
  statistic = c(unname(t_res$statistic), unname(t_student$statistic),
                unname(wilcox_res$statistic)),
  p_value   = c(t_res$p.value, t_student$p.value, wilcox_res$p.value)
)
write.csv(t_test_table, "tables/t2_ttest_results.csv", row.names = FALSE)


# --- 2.ii One-sample proportion test: proportion of male penguins ------------
# H0: proportion of males = 0.50
# H1: proportion of males != 0.50   (two-sided)
cat("\n--- 2.ii One-sample proportion test: proportion of males ---\n")

n_known  <- sum(!is.na(penguins$sex))
n_male   <- sum(penguins$sex == "male",   na.rm = TRUE)
n_female <- sum(penguins$sex == "female", na.rm = TRUE)

cat("\nPenguins with recorded sex:", n_known,
    "(males =", n_male, ", females =", n_female, ")\n")
cat("Sample proportion of males:", round(n_male / n_known, 4), "\n")
cat("Expected males under H0   :", n_known * 0.5, "\n")

# Main test: prop.test (normal approximation with continuity correction)
prop_res <- prop.test(x = n_male, n = n_known, p = 0.5,
                      alternative = "two.sided", correct = TRUE,
                      conf.level = 1 - alpha)
cat("\nMAIN TEST - one-sample proportion test (prop.test):\n")
print(prop_res)

# Robustness check: exact binomial test
cat("\nRobustness check - exact binomial test:\n")
binom_res <- binom.test(x = n_male, n = n_known, p = 0.5)
print(binom_res)

cat("\np-value (prop.test):", fmt_p(prop_res$p.value), "\n")
cat("Decision at alpha =", alpha, ":",
    ifelse(prop_res$p.value < alpha,
           "REJECT H0 - proportion of males differs from 50%.",
           "FAIL TO REJECT H0 - no evidence the proportion differs from 50%."),
    "\n")

prop_table <- data.frame(
  test       = c("prop.test", "binom.test"),
  n_male     = n_male,
  n          = n_known,
  prop_male  = n_male / n_known,
  p_value    = c(prop_res$p.value, binom_res$p.value)
)
write.csv(prop_table, "tables/t2_proportion_results.csv", row.names = FALSE)


# =============================================================================
# 3. TASK 3: CORRELATION ANALYSIS
# =============================================================================
section("TASK 3: CORRELATION ANALYSIS")

# --- 3.1 Correlation matrix (Pearson) ----------------------------------------
cor_mat <- cor(as.data.frame(penguins)[num_vars], use = "complete.obs")
cat("\nPearson correlation matrix (n = 342 complete measurements):\n")
print(round(cor_mat, 3))

write.csv(round(cor_mat, 4), "tables/t3_correlation_matrix.csv")

# Significance test for every pair of variables
pair_list <- combn(num_vars, 2, simplify = FALSE)
cor_tests <- do.call(rbind, lapply(pair_list, function(p) {
  ct <- cor.test(penguins[[p[1]]], penguins[[p[2]]])
  data.frame(
    variable_1 = p[1],
    variable_2 = p[2],
    r          = unname(ct$estimate),
    ci_low     = ct$conf.int[1],
    ci_high    = ct$conf.int[2],
    p_value    = ct$p.value
  )
}))
cor_tests$strength <- cut(abs(cor_tests$r),
                          breaks = c(0, 0.3, 0.5, 0.7, 1),
                          labels = c("weak", "moderate", "fairly strong", "high"),
                          include.lowest = TRUE)
cat("\nCorrelation tests for each pair:\n")
print(transform(cor_tests,
                r = round(r, 3), ci_low = round(ci_low, 3),
                ci_high = round(ci_high, 3),
                p_value = signif(p_value, 3)),
      row.names = FALSE)
write.csv(cor_tests, "tables/t3_correlation_tests.csv", row.names = FALSE)

# Extra check: do the relationships hold WITHIN each species?
within_species_cor <- penguins %>%
  filter(!is.na(bill_length_mm)) %>%
  group_by(species) %>%
  summarise(
    n                    = n(),
    r_billlen_billdepth  = round(cor(bill_length_mm, bill_depth_mm), 3),
    r_flipper_bodymass   = round(cor(flipper_length_mm, body_mass_g), 3),
    .groups = "drop"
  )
cat("\nCorrelations computed separately within each species:\n")
print(as.data.frame(within_species_cor))
write.csv(within_species_cor, "tables/t3_within_species_cor.csv", row.names = FALSE)

# --- 3.2 Heatmap of the correlation matrix -----------------------------------
pretty_names <- c("Bill length", "Bill depth", "Flipper length", "Body mass")
cor_plot_mat <- cor_mat
dimnames(cor_plot_mat) <- list(pretty_names, pretty_names)

cor_long <- as.data.frame(as.table(cor_plot_mat))
names(cor_long) <- c("Var1", "Var2", "r")

p_heat <- ggplot(cor_long, aes(x = Var1, y = Var2, fill = r)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 5) +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
                       midpoint = 0, limits = c(-1, 1),
                       name = "Pearson r") +
  coord_fixed() +
  labs(title = "Correlation matrix of penguin measurements",
       subtitle = "Pearson correlation, n = 342 penguins with complete measurements",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_blank())
print(p_heat)
ggsave("figures/fig3_correlation_heatmap.png", p_heat,
       width = 7.5, height = 6, dpi = 300)

# --- 3.3 Discussion of highly correlated variables ---------------------------
# (See the PDF report for the full written discussion.)
cat("\nPairs with |r| >= 0.70 (highly correlated):\n")
print(cor_tests[abs(cor_tests$r) >= 0.70,
                c("variable_1", "variable_2", "r")], row.names = FALSE)


# =============================================================================
# 4. TASK 4: REGRESSION MODELLING
# =============================================================================
section("TASK 4: REGRESSION MODELLING")

# --- 4.i Simple linear regression: body_mass_g ~ flipper_length_mm -----------
cat("\n--- 4.i Simple linear regression ---\n")
model1 <- lm(body_mass_g ~ flipper_length_mm, data = penguins)
print(summary(model1))

b <- coef(model1)
cat("\nRegression equation:\n")
cat(sprintf("  body_mass_g = %.2f + %.2f * flipper_length_mm\n", b[1], b[2]))
cat("\n95% confidence intervals for the coefficients:\n")
print(round(confint(model1), 3))
cat(sprintf("\nR-squared = %.4f, adjusted R-squared = %.4f\n",
            summary(model1)$r.squared, summary(model1)$adj.r.squared))
cat(sprintf("Residual standard error = %.1f g on %d degrees of freedom\n",
            summary(model1)$sigma, model1$df.residual))

# --- 4.ii Multiple linear regression -----------------------------------------
cat("\n--- 4.ii Multiple linear regression ---\n")
model2 <- lm(body_mass_g ~ flipper_length_mm + bill_length_mm + bill_depth_mm,
             data = penguins)
print(summary(model2))
cat("\n95% confidence intervals for the coefficients:\n")
print(round(confint(model2), 3))

# Multicollinearity: variance inflation factors (VIF = 1 / (1 - R^2_j))
predictors <- c("flipper_length_mm", "bill_length_mm", "bill_depth_mm")
vif_values <- sapply(predictors, function(v) {
  others <- setdiff(predictors, v)
  r2 <- summary(lm(reformulate(others, response = v), data = penguins))$r.squared
  1 / (1 - r2)
})
cat("\nVariance inflation factors (values above 5 would be a concern):\n")
print(round(vif_values, 2))

# Model comparison
cat("\nF-test comparing the simple model with the multiple model:\n")
model_anova <- anova(model1, model2)
print(model_anova)

model_comparison <- data.frame(
  model         = c("Simple (flipper length)",
                    "Multiple (flipper + bill length + bill depth)"),
  predictors    = c(1, 3),
  r_squared     = c(summary(model1)$r.squared,     summary(model2)$r.squared),
  adj_r_squared = c(summary(model1)$adj.r.squared, summary(model2)$adj.r.squared),
  rse_grams     = c(summary(model1)$sigma,         summary(model2)$sigma),
  AIC           = c(AIC(model1), AIC(model2)),
  BIC           = c(BIC(model1), BIC(model2))
)
cat("\nModel comparison:\n")
print(transform(model_comparison,
                r_squared = round(r_squared, 4),
                adj_r_squared = round(adj_r_squared, 4),
                rse_grams = round(rse_grams, 1),
                AIC = round(AIC, 1), BIC = round(BIC, 1)),
      row.names = FALSE)
write.csv(model_comparison, "tables/t4_model_comparison.csv", row.names = FALSE)

coef_table <- function(m) {
  ct <- as.data.frame(summary(m)$coefficients)
  ci <- as.data.frame(confint(m))
  out <- data.frame(term = rownames(ct), ct, ci, row.names = NULL)
  names(out) <- c("term", "estimate", "std_error", "t_value", "p_value",
                  "ci_low", "ci_high")
  out
}
write.csv(coef_table(model1), "tables/t4_model1_coefficients.csv", row.names = FALSE)
write.csv(coef_table(model2), "tables/t4_model2_coefficients.csv", row.names = FALSE)

# --- 4.iii Diagnostic checks for the multiple regression ---------------------
mf <- model.frame(model2)

# 1. Normality of the residuals (Shapiro-Wilk)
cat("\n1. Shapiro-Wilk test on the residuals of the multiple regression:\n")
shapiro_resid <- shapiro.test(residuals(model2))
print(shapiro_resid)

# 2. Constant spread: Breusch-Pagan (studentised) test, calculated by hand
aux <- lm(residuals(model2)^2 ~ flipper_length_mm + bill_length_mm +
            bill_depth_mm, data = mf)
bp_stat <- nobs(model2) * summary(aux)$r.squared
bp_p    <- pchisq(bp_stat, df = 3, lower.tail = FALSE)
cat(sprintf("\n2. Breusch-Pagan test: statistic = %.3f, p-value = %s\n",
            bp_stat, fmt_p(bp_p)))

# 3. Influence: Cook's distance and large standardised residuals
max_cook   <- max(cooks.distance(model2))
n_outliers <- sum(abs(rstandard(model2)) > 3)
cat(sprintf("\n3. Largest Cook's distance = %.3f (concern if > 0.5)\n", max_cook))
cat("   Penguins with |standardised residual| > 3:", n_outliers, "\n")

# 4. Linearity: RESET test (does adding the squared fitted values help?)
mf$fit2   <- fitted(model2)^2
reset_fit <- lm(body_mass_g ~ flipper_length_mm + bill_length_mm +
                  bill_depth_mm + fit2, data = mf)
reset_p   <- anova(model2, reset_fit)$"Pr(>F)"[2]
cat(sprintf("\n4. RESET test for non-linearity: p-value = %s\n", fmt_p(reset_p)))

png("figures/fig4_diagnostics_model2.png", width = 1800, height = 1800, res = 200)
par(mfrow = c(2, 2))
plot(model2)
par(mfrow = c(1, 1))
dev.off()

# --- 4.iv EXTRA (not required): does species change the picture? -------------
# Adding species shows how much of the flipper-body mass relationship is
# explained by differences BETWEEN species. Used in the limitations section.
model_species <- lm(body_mass_g ~ flipper_length_mm + species, data = penguins)
cat("\nExtra model: body_mass_g ~ flipper_length_mm + species\n")
print(summary(model_species))


# =============================================================================
# 5. TASK 5: VISUALIZATION
# =============================================================================
section("TASK 5: VISUALIZATION")

pen_fl_bm <- penguins[!is.na(penguins$flipper_length_mm) &
                        !is.na(penguins$body_mass_g), ]

# --- 5.1 Scatter plot of body mass vs flipper length with regression line ----
eq_label <- sprintf("Body mass = %.1f + %.2f x Flipper length\nR-squared = %.3f",
                    b[1], b[2], summary(model1)$r.squared)

p_scatter <- ggplot(pen_fl_bm, aes(x = flipper_length_mm, y = body_mass_g)) +
  geom_point(aes(colour = species), alpha = 0.75, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "black", linewidth = 0.9) +
  annotate("text", x = min(pen_fl_bm$flipper_length_mm) + 1,
           y = max(pen_fl_bm$body_mass_g), label = eq_label,
           hjust = 0, vjust = 1, size = 4) +
  scale_colour_manual(values = species_cols, name = "Species") +
  labs(title = "Body mass vs flipper length",
       subtitle = "Black line = least-squares regression line (with 95% confidence band)",
       x = "Flipper length (mm)", y = "Body mass (g)") +
  theme_minimal(base_size = 12)
print(p_scatter)
ggsave("figures/fig1_scatter_flipper_body_mass.png", p_scatter,
       width = 8, height = 6, dpi = 300)

# --- 5.2 Boxplot of body mass across species ---------------------------------
pen_bm_sp <- penguins[!is.na(penguins$body_mass_g), ]

p_box <- ggplot(pen_bm_sp, aes(x = species, y = body_mass_g, fill = species)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, alpha = 0.35, size = 1.3) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = species_cols) +
  guides(fill = "none") +
  labs(title = "Distribution of body mass by species",
       subtitle = "Diamond = group mean; box = median and interquartile range",
       x = "Species", y = "Body mass (g)") +
  theme_minimal(base_size = 12)
print(p_box)
ggsave("figures/fig2_boxplot_body_mass_species.png", p_box,
       width = 8, height = 6, dpi = 300)

# --- 5.3 Histogram of body mass with normal curve overlay --------------------
bm <- penguins$body_mass_g[!is.na(penguins$body_mass_g)]

p_hist <- ggplot(data.frame(body_mass_g = bm), aes(x = body_mass_g)) +
  geom_histogram(aes(y = after_stat(density)), bins = 25,
                 fill = "#8DA0CB", colour = "white") +
  stat_function(fun = dnorm, args = list(mean = mean(bm), sd = sd(bm)),
                colour = "#B2182B", linewidth = 1) +
  labs(title = "Histogram of body mass with fitted normal curve",
       subtitle = sprintf("Normal curve uses mean = %.0f g and SD = %.0f g",
                          mean(bm), sd(bm)),
       x = "Body mass (g)", y = "Density") +
  theme_minimal(base_size = 12)
print(p_hist)
ggsave("figures/fig5_histogram_body_mass.png", p_hist,
       width = 8, height = 6, dpi = 300)

cat("\nSkewness check for body mass (positive = longer right tail):",
    round(mean((bm - mean(bm))^3) / sd(bm)^3, 3), "\n")
cat("Shapiro-Wilk test on body mass:\n")
print(shapiro.test(bm))


# =============================================================================
# 6. SESSION INFORMATION (for reproducibility)
# =============================================================================
section("SESSION INFORMATION")
print(sessionInfo())

cat("\nDone. Figures are in 'figures/' and tables are in 'tables/'.\n")