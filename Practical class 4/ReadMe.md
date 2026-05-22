# Practice 4: Correlation and Regression Analysis of Clinical Hormone Data

> Spearman Correlation · Permutation Testing · Regression Model Comparison · Logistic Regression · ROC Analysis

---

## Project Overview

This project applies a structured correlation and regression analysis pipeline to a clinical dataset. Starting from the cleaned, imputed dataset produced in Practice 1, the work progresses through Spearman correlation with permutation-based significance testing, multiple regression model fitting and BIC-based comparison, and binomial logistic regression to predict tumour presence using hormone variables — all implemented in R from a single script.

---

## Data and R Environment

### Dataset

```r
data <- read.csv("data_for_analysis.csv")
data$outcome <- as.factor(data$outcome)
```

**Source:** Cleaned and imputed output from Practice 1.

| Property | Detail |
|---|---|
| Object name | `data` |
| Records | 1,148 rows |
| Variables | 31 columns |
| Outcome = 0 (no tumour) | 987 records |
| Outcome = 1 (tumour) | 160 records |
| Key variable groups | `hormone1`–`hormone10_generated`, `lipids1`–`4`, `antioxidant1`–`5`, `lipid_pero1`–`5`, `carb_metabolism`, `outcome` |

### R Environment

| Component | Detail |
|---|---|
| IDE | RStudio (Windows) |
| R version | 4.3.3 |
| `wPerm` (v1.0.1) | Permutation-based correlation significance via `perm.relation()` — installed via `remotes::install_version()` |
| `pROC` | ROC curve and AUC computation via `roc()` and `auc()` |
| `ggplot2` | Correlation heatmap visualisation |

> **Note:** `wPerm` is not available on CRAN for newer R versions. Install with:
> ```r
> install.packages("remotes")
> remotes::install_version("wPerm", version = "1.0.1")
> ```

---

## Initial Exploration — Normality Testing and Lipids Correlation (Worked Example)

Before the main tasks, the script verifies the normality of `lipids1` and `lipids2` and runs a worked Spearman correlation example as a pipeline check:

```r
shapiro.test(data$lipids1)
shapiro.test(data$lipids2)

hist(data$lipids1)
qqnorm(data$lipids1)

spearman_result <- cor.test(data$lipids1, data$lipids2, method = "spearman")
print(spearman_result)
```

A permutation Spearman test (10,000 permutations) is then run for `lipids1` against `lipids2`, `lipids3`, and `lipids4`:

```r
target_vars <- c("lipids2", "lipids3", "lipids4")

for (var in target_vars) {
  perm_spearman <- perm.relation(
    x      = data$lipids1,
    y      = data[[var]],
    method = "spearman",
    R      = 10000
  )
  results <- rbind(results, data.frame(
    variable      = var,
    spearman_corr = perm_spearman$Observed,
    s_p_value     = perm_spearman$p.value
  ))
}
print(results)
```

The hormone variable vector used throughout all subsequent tasks:

```r
hormone_cols <- c("hormone1", "hormone2", "hormone3", "hormone4",
                  "hormone5", "hormone6", "hormone7", "hormone8",
                  "hormone10_generated")
```

---

## Task 1 — Spearman Correlation Matrix (All Hormones)

Since all hormone variables were confirmed as non-normally distributed in Practice 3 (Shapiro–Wilk, p < 0.05), Spearman rank correlation is used throughout. A full 9×9 Spearman correlation matrix is computed across all hormone pairs:

```r
cor_matrix <- cor(data[, hormone_cols],
                  use    = "pairwise.complete.obs",
                  method = "spearman")
print(cor_matrix)
```

**Finding:** The strongest correlation was between `hormone3` and `hormone4` (rho = 0.584), indicating a moderate positive relationship. Most other hormone pairs showed weak or negligible correlations.

---

## Task 2 — Permutation-Based Significance Testing (All Hormone Pairs)

A permutation test (1,000 permutations, `R = 1000`) is applied to all 36 unique hormone pairs using `perm.relation()` from `wPerm`. Results include the observed Spearman rho and a permutation-based p-value:

```r
hormone_results <- data.frame()

for (i in 1:(length(hormone_cols) - 1)) {
  for (j in (i + 1):length(hormone_cols)) {
    perm <- perm.relation(
      x      = data[[hormone_cols[i]]],
      y      = data[[hormone_cols[j]]],
      method = "spearman",
      R      = 1000
    )
    hormone_results <- rbind(hormone_results, data.frame(
      Variable1    = hormone_cols[i],
      Variable2    = hormone_cols[j],
      Spearman_rho = round(perm$Observed, 5),
      p_value      = perm$p.value
    ))
  }
}

write.csv(hormone_results, "task1_2_hormone_correlation_permutation.csv", row.names = FALSE)
```

**Significant pairs (p < 0.05):**

| Variable 1 | Variable 2 | Spearman rho | Direction |
|---|---|---|---|
| hormone3 | hormone4 | 0.584 | Positive (strongest) |
| hormone4 | hormone7 | −0.325 | Negative |
| hormone5 | hormone8 | 0.259 | Positive |
| hormone3 | hormone5 | 0.241 | Positive |
| hormone7 | hormone10_generated | 0.200 | Positive |

| Output | File |
|---|---|
| Hormone pairwise Spearman + permutation p-values | `task1_2_hormone_correlation_permutation.csv` |

---

## Task 3 — Regression Analysis Between Hormone Variables

The pair with the strongest Spearman correlation — `hormone3` (predictor) and `hormone4` (outcome), rho = 0.584, p = 0.002 — is selected for regression analysis. The data is ordered by `hormone3` before fitting:

```r
# Selected pair: hormone3 (predictor) and hormone4 (outcome)
# Reason: strongest Spearman correlation among all hormone pairs (rho = 0.584, p = 0.002)

df_h <- data
df_h <- df_h[order(df_h$hormone3), ]

# Linear regression
model_h_linear <- lm(hormone4 ~ hormone3, data = df_h)

# 2nd-degree polynomial
model_h_2 <- lm(hormone4 ~ poly(hormone3, 2), data = df_h)

# 3rd-degree polynomial
model_h_3 <- lm(hormone4 ~ poly(hormone3, 3), data = df_h)

# Exponential dependence (log-transformed response)
model_h_exp <- lm(log(hormone4) ~ hormone3, data = df_h)

# Logarithmic dependence (exp-transformed response)
model_h_log <- lm(exp(hormone4) ~ hormone3, data = df_h)
```

**Finding:** The exponential model (`log(hormone4) ~ hormone3`) achieved the highest R-squared (0.297). The logarithmic model (`exp(hormone4) ~ hormone3`) performed very poorly due to the extreme values produced by `exp(hormone4)`.

---

## Task 4 — Model Selection via BIC

BIC is computed for all five regression models and ranked in ascending order (lower = better fit):

```r
rezult_h <- data.frame(
  model = c("model_h_linear", "model_h_2", "model_h_3",
            "model_h_exp",    "model_h_log"),
  BIC_value = c(BIC(model_h_linear), BIC(model_h_2), BIC(model_h_3),
                BIC(model_h_exp),    BIC(model_h_log))
)

rezult_h <- rezult_h[order(rezult_h$BIC_value), ]
print(rezult_h)
write.csv(rezult_h, "task4_hormone_BIC_comparison.csv", row.names = FALSE)
```

**BIC Results (ranked):**

| Rank | Model | Formula | BIC Value |
|---|---|---|---|
| 1 ✓ | `model_h_exp` | `log(hormone4) ~ hormone3` | 2,049.49 |
| 2 | `model_h_2` | `hormone4 ~ poly(hormone3, 2)` | 8,476.35 |
| 3 | `model_h_3` | `hormone4 ~ poly(hormone3, 3)` | 8,477.46 |
| 4 | `model_h_linear` | `hormone4 ~ hormone3` | 8,525.91 |
| 5 | `model_h_log` | `exp(hormone4) ~ hormone3` | 296,869.92 |

**Conclusion:** The exponential model (`log(hormone4) ~ hormone3`) is the best-fitting model, confirmed by both the lowest BIC (2,049.49) and the highest R-squared (0.297).

| Output | File |
|---|---|
| BIC ranking of all 5 hormone regression models | `task4_hormone_BIC_comparison.csv` |

---

## Task 5 — Logistic Regression to Predict Binary Outcome

Three binomial logistic regression models are fitted to predict `outcome` (0 = no tumour, 1 = tumour) using hormone variables:

```r
# Simple model: one predictor
model_logit_h1 <- glm(outcome ~ hormone1,
                      data = data, family = binomial)

# Two predictors: strongest correlated pair
model_logit_h2 <- glm(outcome ~ hormone3 + hormone4,
                      data = data, family = binomial)

# Full model: all nine hormones
model_logit_hall <- glm(outcome ~ hormone1 + hormone2 + hormone3 +
                          hormone4 + hormone5 + hormone6 + hormone7 +
                          hormone8 + hormone10_generated,
                        data = data, family = binomial)
```

### AIC/BIC Model Comparison

```r
aic_bic <- data.frame(
  Model = c("model_logit_h1", "model_logit_h2", "model_logit_hall"),
  AIC   = c(AIC(model_logit_h1), AIC(model_logit_h2), AIC(model_logit_hall)),
  BIC   = c(BIC(model_logit_h1), BIC(model_logit_h2), BIC(model_logit_hall))
)
print(aic_bic)
write.csv(aic_bic, "task5_logistic_AIC_BIC.csv", row.names = FALSE)
```

| Model | Predictors | AIC | BIC |
|---|---|---|---|
| `model_logit_h1` | `hormone1` | 928.83 | 938.92 |
| `model_logit_h2` | `hormone3 + hormone4` | 929.69 | 944.82 |
| `model_logit_hall` ✓ | All 9 hormones | **919.31** | 969.76 |

**Best by AIC:** `model_logit_hall` (AIC = 919.31).

### Stepwise Variable Selection

```r
step_model_h <- step(model_logit_hall, direction = "both")
summary(step_model_h)
```

Stepwise selection (AIC, `direction = "both"`) reduced the full model to four predictors: `hormone1`, `hormone2`, `hormone5`, and `hormone8`, with AIC = 913.12.

**`hormone8` was the strongest and only statistically significant predictor** (p = 0.00122), with OR = 0.996 (95% CI: 0.994–0.998), indicating that higher `hormone8` levels are associated with lower odds of tumour presence.

### Probability Prediction and Confusion Matrix

```r
data$pred_prob_h  <- predict(model_logit_hall, type = "response")
data$pred_class_h <- ifelse(data$pred_prob_h > 0.5, 1, 0)

conf_matrix <- as.data.frame(table(Actual    = data$outcome,
                                   Predicted = data$pred_class_h))
write.csv(conf_matrix, "task5_confusion_matrix.csv", row.names = FALSE)
```

**Finding:** The model predicted all cases as non-tumour (class 0) at the 0.5 threshold, reflecting the class imbalance in the dataset (987 non-tumour vs 160 tumour cases). This is a known limitation of logistic regression applied to heavily imbalanced binary outcomes without resampling or threshold adjustment.

### Odds Ratios and 95% Confidence Intervals

```r
odds_ratios <- as.data.frame(
  exp(cbind(OR = coef(model_logit_hall), confint(model_logit_hall)))
)
write.csv(odds_ratios, "task5_odds_ratios.csv", row.names = TRUE)
```

### ROC Curve and AUC

```r
library(pROC)
roc_curve_h <- roc(data$outcome, data$pred_prob_h)
plot(roc_curve_h, main = "ROC Curve - Hormone Logistic Regression")
auc(roc_curve_h)
```

**AUC = 0.6262** — moderate discriminative ability above random chance (AUC = 0.5). The model distinguishes tumour from non-tumour cases better than chance, but not strongly, due to class imbalance and limited predictive power of hormone variables alone.

| Output | File |
|---|---|
| AIC/BIC for all 3 hormone logistic models | `task5_logistic_AIC_BIC.csv` |
| Confusion matrix — full hormone model | `task5_confusion_matrix.csv` |
| Odds ratios + 95% CI — full hormone model | `task5_odds_ratios.csv` |
| ROC curve plot | `09_roc_curve_hormone_logistic.png` |

---

## Supplementary — Lipids Regression and Logistic Analysis

The script also performs a parallel regression pipeline on `lipids1` (outcome) and `lipids2` (predictor) for comparative reference, along with lipids-based logistic regression evaluated using a ROC curve:

```r
df <- data[order(data$lipids1), ]

model_linear <- lm(lipids1 ~ lipids2, data = df)
model_2      <- lm(lipids1 ~ poly(lipids2, 2), data = df)
model_3      <- lm(lipids1 ~ poly(lipids2, 3), data = df)
model_exp    <- lm(log(lipids1) ~ lipids2, data = df)
model_log    <- lm(exp(lipids1) ~ lipids2, data = df)

rezult <- data.frame(
  model     = c("model_linear", "model_2", "model_3", "model_exp", "model_log"),
  BIC_value = c(BIC(model_linear), BIC(model_2), BIC(model_3),
                BIC(model_exp),    BIC(model_log))
)
rezult <- rezult[order(rezult$BIC_value), ]
```

Scatter plot with OLS regression line and fitted model overlay:

```r
data <- data[order(data$lipids2), ]
plot(data$lipids2, data$lipids1)
lines(data$lipids2, data$lipids1, col = "blue")
abline(lm(data$lipids1 ~ data$lipids2), col = "red")

plot(df$lipids2, df$lipids1)
lines(df$lipids2, fitted(model_linear), col = "blue")
```

Lipids logistic regression and ROC evaluation:

```r
model_logit_1   <- glm(outcome ~ lipids1, data = data, family = binomial)
model_logit_2   <- glm(outcome ~ lipids1 + lipids2, data = data, family = binomial)
model_logit_all <- glm(outcome ~ lipids1 + lipids2 + lipids3 + lipids4,
                       data = data, family = binomial)

data$pred_prob  <- predict(model_logit_2, type = "response")
data$pred_class <- ifelse(data$pred_prob > 0.5, 1, 0)

roc_curve <- roc(data$outcome, data$pred_prob)
plot(roc_curve, main = "ROC-Curve")
auc(roc_curve)

step_model <- step(model_logit_all, direction = "both")
exp(cbind(OR = coef(model_logit_2), confint(model_logit_2)))
```

---

## Complete Output Summary

### Plots (`practice4_plots/`)

| # | File | Description |
|---|---|---|
| 01 | `01_histogram_lipids1.png` | Normality check — histogram of `lipids1` |
| 02 | `02_qqplot_lipids1.png` | Normality check — Q-Q plot of `lipids1` |
| 03 | `03_scatter_lipids2_vs_lipids1_with_regression.png` | Scatter plot with OLS regression line |
| 04 | `04_linear_regression_lipids1_vs_lipids2.png` | Linear model fitted line |
| 05 | `05_all_regression_models_lipids1_vs_lipids2.png` | Linear, Poly-2, Poly-3 overlay |
| 06 | `06_roc_curve_lipids_logistic.png` | ROC curve — lipids logistic model |
| 07 | `07_scatter_hormone3_vs_hormone4.png` | `hormone3` vs `hormone4` scatter with linear fit |
| 08 | `08_all_regression_models_hormone3_vs_hormone4.png` | All 5 regression models overlay — hormones |
| 09 | `09_roc_curve_hormone_logistic.png` | ROC curve — full hormone logistic model (AUC = 0.6262) |
| 10 | `10_BIC_comparison_lipids_models.png` | BIC bar chart — lipids regression models |
| 11 | `11_BIC_comparison_hormone_models.png` | BIC bar chart — hormone regression models |
| 12 | `12_spearman_heatmap_all_hormones.png` | Spearman correlation heatmap — all 9 hormones |

### CSV Files (`practice4_csv/`)

| File | Contents |
|---|---|
| `lipids1_permutation_spearman.csv` | Permutation Spearman: `lipids1` vs `lipids2`, `3`, `4` (R = 10,000) |
| `task1_2_hormone_correlation_permutation.csv` | All 36 hormone pairwise Spearman rho + permutation p-values (R = 1,000) |
| `task4_hormone_BIC_comparison.csv` | BIC ranking of 5 hormone regression models |
| `lipids_BIC_comparison.csv` | BIC ranking of 5 lipid regression models |
| `task5_logistic_AIC_BIC.csv` | AIC and BIC for all 3 hormone logistic models |
| `task5_confusion_matrix.csv` | Confusion matrix — full hormone logistic model at threshold 0.5 |
| `task5_odds_ratios.csv` | OR + 95% CI for all predictors in `model_logit_hall` |
