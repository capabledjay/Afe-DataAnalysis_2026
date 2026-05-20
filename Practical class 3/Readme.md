
## Project Overview

This project applies a structured statistical analysis pipeline to a clinical dataset. Starting from a cleaned, imputed dataset produced in Practice 1, the work progresses through descriptive statistics, normality and variance testing, multi-test hypothesis comparison, and exploratory correlation analysis — all implemented in R and reproducible from a single script.

---

## Data and R Environment

### Dataset

```r
data_for_analysis <- read.csv("data_for_analysis.csv")
```

**Source:** Cleaned and imputed output from Practice 1.

| Property | Detail |
|---|---|
| Object name | `data_for_analysis` |
| Records | 1,148 rows |
| Variables | 31 columns |
| Outcome = 0 (no tumour) | 987 records |
| Outcome = 1 (tumour) | 160 records |
| Key variable groups | `hormone1`–`hormone10_generated`, `lipids1`–`5`, `antioxidant1`–`5`, `lipid_pero1`–`5`, `carb_metabolism`, `outcome`, factor variables |

### R Environment

| Component | Detail |
|---|---|
| IDE | RStudio (Windows) |
| R version | 4.3.3 |
| `gtsummary` / `cardx` | Publication-ready descriptive statistics tables via `tbl_summary()` |
| `car` | Levene's Test for homogeneity of variance via `car::leveneTest()` |
| `lawstat` | Brunner–Munzel test via `brunner.munzel.test()` |
| `ggplot2` + `reshape2` | Spearman correlation heatmap construction and data reshaping via `melt()` |
| `DataExplorer` | Automated EDA HTML report via `create_report()` |

---

## Initial Exploration — lipids1 (Worked Example)

Before the formal task loop, the script performs a worked example on `lipids1` to verify the testing workflow:

```r
value_outcome1 <- data_for_analysis[data_for_analysis$outcome == "1", ]$lipids1
hist(value_outcome1, col = "lightblue")
qqnorm(value_outcome1, main = "Q-Q Plot")
qqline(value_outcome1, col = "red", lwd = 2)
shapiro.test(value_outcome1)

value_outcome0 <- data_for_analysis[data_for_analysis$outcome == "0", ]$lipids1
hist(value_outcome0, col = "lightgreen")
qqnorm(value_outcome0, main = "Q-Q Plot")
qqline(value_outcome0, col = "red", lwd = 2)
shapiro.test(value_outcome0)
```

A Levene's test and Brunner–Munzel test, alongside `t.test()` and `wilcox.test()`, are run on the two `lipids1` groups as a reference before the hormone-level loop begins:

```r
car::leveneTest(lipids1 ~ outcome, data = data_for_analysis)
brunner.munzel.test(group1, group2)
t.test(group1, group2)
wilcox.test(group1, group2)
```

The hormone variable vector used throughout all subsequent tasks:

```r
hormone_cols <- c("hormone1", "hormone2", "hormone3", "hormone4",
                  "hormone5", "hormone6", "hormone7", "hormone8",
                  "hormone10_generated")
```

---

## Task 1 — Descriptive Statistics Table by Group (All Hormones)

A publication-ready descriptive statistics table is created for all nine hormone variables grouped by `outcome` using `tbl_summary()` from `gtsummary`. Because all hormone variables were confirmed as non-normally distributed (see Task 2), parameters are reported as **median and IQR**. A Wilcoxon-based p-value is appended via `add_p()`.

```r
tbl_summary(
  data_for_analysis[, c("outcome", hormone_cols)],
  by = outcome,
  statistic = list(all_continuous() ~ "{median} ({p25}, {p75})"),
  digits = all_continuous() ~ 3
) %>%
  add_p() %>%
  bold_labels() %>%
  modify_caption("**Table 3. Descriptive Statistics of Hormones by Outcome Group with p-values**") %>%
  as_tibble() %>%
  write.csv("task1_hormone_descriptive_stats_table.csv", row.names = FALSE)
```

| Item | Detail |
|---|---|
| Statistic reported | Median (Q25, Q75) |
| Grouping variable | `outcome` (0 = no tumour, 1 = tumour) |
| P-value method | Wilcoxon rank-sum via `add_p()` |
| Output file | `task1_hormone_descriptive_stats_table.csv` |

---

## Task 2 — Shapiro–Wilk Test and Levene's Test (All Hormones)

### Shapiro–Wilk Normality Test

Run for all nine hormone variables in both outcome groups using a nested for-loop. Collects W-statistic, p-value, and a binary normality flag per hormone per group:

```r
sw_results <- data.frame()
for (h in hormone_cols) {
  for (grp in c("0", "1")) {
    vals <- data_for_analysis[[h]][data_for_analysis$outcome == grp]
    vals <- vals[!is.na(vals)]
    sw   <- shapiro.test(vals)
    sw_results <- rbind(sw_results, data.frame(
      Hormone = h, Group = grp,
      W       = round(sw$statistic, 5),
      p_value = round(sw$p.value, 5),
      Normal  = ifelse(sw$p.value >= 0.05, "Yes", "No")))
  }
}
write.csv(sw_results, "shapiro_wilk_test_results.csv", row.names = FALSE)
```

> **Finding:** All hormone variables were non-normally distributed (p < 0.05) in both outcome groups. This justifies the use of non-parametric methods throughout Tasks 1, 4, and 5.

### Levene's Test for Homogeneity of Variance

Run for all nine hormone variables with `center = median`, in line with non-normal data:

```r
lev_results <- data.frame()
for (h in hormone_cols) {
  lev <- car::leveneTest(data_for_analysis[[h]] ~ data_for_analysis$outcome,
                         center = median)
  lev_results <- rbind(lev_results, data.frame(
    Hormone   = h,
    F_value   = round(lev$`F value`[1], 4),
    p_value   = round(lev$`Pr(>F)`[1], 5),
    Equal_Var = ifelse(lev$`Pr(>F)`[1] >= 0.05, "Yes", "No")))
}
write.csv(lev_results, "levene_test_results.csv", row.names = FALSE)
```

> **Finding:** All hormone variables showed equal variances between outcome groups (p > 0.05). This finding, combined with non-normality, directly informs test selection in Task 4.

| Output | File |
|---|---|
| Shapiro–Wilk results | `shapiro_wilk_test_results.csv` |
| Levene's test results | `levene_test_results.csv` |

---

## Task 3 — Histograms and Q-Q Plots (All Hormones)

### Histograms

Side-by-side histograms (`par(mfrow = c(1, 2))`) for each hormone, comparing outcome 0 (`lightblue`) vs outcome 1 (`salmon`), with 30 breaks:

```r
for (h in hormone_cols) {
  par(mfrow = c(1, 2))
  g0 <- data_for_analysis[[h]][data_for_analysis$outcome == "0"]
  g1 <- data_for_analysis[[h]][data_for_analysis$outcome == "1"]
  hist(g0, col = "lightblue", main = paste(h, "| Outcome = 0"), xlab = h, breaks = 30)
  hist(g1, col = "salmon",    main = paste(h, "| Outcome = 1"), xlab = h, breaks = 30)
}
par(mfrow = c(1, 1))
```

### Q-Q Plots

Outcome 0 points in `steelblue`, outcome 1 in `darkorange`, both with a red reference line (`lwd = 2`):

```r
for (h in hormone_cols) {
  par(mfrow = c(1, 2))
  qqnorm(g0, main = paste("Q-Q:", h, "| Outcome = 0"), col = "steelblue")
  qqline(g0, col = "red", lwd = 2)
  qqnorm(g1, main = paste("Q-Q:", h, "| Outcome = 1"), col = "darkorange")
  qqline(g1, col = "red", lwd = 2)
}
par(mfrow = c(1, 1))
```

> **Finding:** All hormone distributions deviated clearly from the Q-Q reference line, confirming non-normality. Most hormones were right-skewed. `hormone10_generated` showed the most extreme deviation in outcome group 0, with a near-zero concentration of values.

| Plot type | Output files |
|---|---|
| Histograms | `05_histogram_hormone1.png` → `13_histogram_hormone10_generated.png` |
| Q-Q plots | `14_qqplot_hormone1.png` → `22_qqplot_hormone10_generated.png` |

---

## Task 4 — Hypothesis Testing: Brunner–Munzel, t-test, and Wilcoxon (All Hormones)

Three statistical tests are applied to each hormone variable comparing outcome group 0 vs outcome group 1. NA values are removed before testing. Results are compiled alongside median (IQR) for each group:

```r
test_results <- data.frame()
for (h in hormone_cols) {
  g0 <- data_for_analysis[[h]][data_for_analysis$outcome == "0"]
  g1 <- data_for_analysis[[h]][data_for_analysis$outcome == "1"]
  g0 <- g0[!is.na(g0)]; g1 <- g1[!is.na(g1)]

  bm <- brunner.munzel.test(g0, g1)
  tt <- t.test(g0, g1)
  wt <- wilcox.test(g0, g1)

  test_results <- rbind(test_results, data.frame(
    Hormone         = h,
    Median_IQR_G0   = paste0(round(median(g0), 3), " (", round(quantile(g0, 0.25), 3), ", ", round(quantile(g0, 0.75), 3), ")"),
    Median_IQR_G1   = paste0(round(median(g1), 3), " (", round(quantile(g1, 0.25), 3), ", ", round(quantile(g1, 0.75), 3), ")"),
    p_BrunnerMunzel = round(bm$p.value, 5),
    p_ttest         = round(tt$p.value, 5),
    p_Wilcoxon      = round(wt$p.value, 5)))
}
write.csv(test_results, "task4_descriptive_with_pvalues.csv", row.names = FALSE)
```

### Test Selection Conclusion

```r
# The Wilcoxon rank-sum test is applicable.
# Because the data is non-normal (Shapiro-Wilk, p < 0.05)
# and variances are equal between groups (Levene's test, p > 0.05).
```

| Test | Applicable? | Reason |
|---|---|---|
| Brunner–Munzel | Valid | Robust to unequal distributions; no normality assumption |
| Welch t-test | Not preferred | Assumes normality — violated here |
| **Wilcoxon rank-sum ✓** | **Recommended** | Non-normal data + equal variances between groups |

| Output | File |
|---|---|
| All test p-values + medians | `task4_descriptive_with_pvalues.csv` |

---

## Task 5 — Spearman Correlation Heatmaps (All Hormones by Group)

Spearman correlation is used throughout because all hormone variables are non-normally distributed (confirmed in Task 2). Separate correlation matrices are computed for each outcome group using pairwise complete observations, then visualised as annotated heatmaps.

### Correlation Matrices

```r
df_g0  <- data_for_analysis[data_for_analysis$outcome == "0", hormone_cols]
cor_g0 <- cor(df_g0, use = "pairwise.complete.obs", method = "spearman")

df_g1  <- data_for_analysis[data_for_analysis$outcome == "1", hormone_cols]
cor_g1 <- cor(df_g1, use = "pairwise.complete.obs", method = "spearman")
```

### Heatmap Construction

Each matrix is reshaped with `melt()` then plotted with `geom_tile()` and `geom_text()`. Colour scale: `steelblue` (−1) → white (0) → `firebrick` (+1):

```r
melted_g0 <- melt(cor_g0)
p_g0 <- ggplot(melted_g0, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 2)), size = 3) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman r") +
  labs(title = "Hormone Correlation Heatmap - Outcome 0 (No Tumour)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

### Saving Outputs

```r
ggsave("task5_correlation_heatmap_outcome0.png", plot = p_g0, width = 8, height = 6, dpi = 300)
ggsave("task5_correlation_heatmap_outcome1.png", plot = p_g1, width = 8, height = 6, dpi = 300)
write.csv(cor_g0, "task5_correlation_matrix_outcome0.csv", row.names = TRUE)
write.csv(cor_g1, "task5_correlation_matrix_outcome1.csv", row.names = TRUE)
```

| Output type | File |
|---|---|
| Heatmap — outcome 0 | `23_task5_correlation_heatmap_outcome0.png` |
| Heatmap — outcome 1 | `24_task5_correlation_heatmap_outcome1.png` |
| Correlation matrix — outcome 0 | `task5_correlation_matrix_outcome0.csv` |
| Correlation matrix — outcome 1 | `task5_correlation_matrix_outcome1.csv` |

---

## Automated EDA Report

```r
library(DataExplorer)
create_report(
  data        = data_for_analysis,
  output_file = "EDA_Report.html",
  output_dir  = getwd(),
  report_title = "EDA Report"
)
```

| Output | File |
|---|---|
| Full EDA HTML report | `EDA_Report.html` (working directory) |

---

## Complete Output Summary

### Plots (`practice3_plots/`)

| # | File | Description |
|---|---|---|
| 01 | `01_histogram_lipids1_outcome1.png` | lipids1 histogram, outcome = 1 |
| 02 | `02_qqplot_lipids1_outcome1.png` | lipids1 Q-Q plot, outcome = 1 |
| 03 | `03_histogram_lipids1_outcome0.png` | lipids1 histogram, outcome = 0 |
| 04 | `04_qqplot_lipids1_outcome0.png` | lipids1 Q-Q plot, outcome = 0 |
| 05–13 | `05–13_histogram_hormone[n].png` | Histograms for all 9 hormones (side-by-side G0/G1) |
| 14–22 | `14–22_qqplot_hormone[n].png` | Q-Q plots for all 9 hormones (side-by-side G0/G1) |
| 23 | `23_task5_correlation_heatmap_outcome0.png` | Spearman heatmap — no tumour group |
| 24 | `24_task5_correlation_heatmap_outcome1.png` | Spearman heatmap — tumour group |

### CSV Files (`practice3_csv/`)

| File | Contents |
|---|---|
| `task1_hormone_descriptive_stats_table.csv` | Median (IQR) + Wilcoxon p-value per hormone by group |
| `shapiro_wilk_test_results.csv` | W, p-value, normality flag per hormone per group |
| `levene_test_results.csv` | F-value, p-value, equal-variance flag per hormone |
| `task4_descriptive_with_pvalues.csv` | Median (IQR) + BM / t-test / Wilcox p-values per hormone |
| `task5_correlation_matrix_outcome0.csv` | 9×9 Spearman r matrix — no tumour group |
| `task5_correlation_matrix_outcome1.csv` | 9×9 Spearman r matrix — tumour group |
