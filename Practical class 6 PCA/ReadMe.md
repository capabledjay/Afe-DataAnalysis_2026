# Practice 6: PCA and Correlation Analysis of Plant Morphometric Data

> Data Standardization · Spearman Correlation · Principal Component Analysis · 3D Visualization

---

## Project Overview

This project applies multivariate dimensionality reduction to a plant morphometric dataset collected across 6 sites. The work progresses through min-max standardization, pairwise Spearman correlation analysis, and Principal Component Analysis (PCA) with 2D and 3D biplot visualization to explore morphological variation and group separation.

---

## Data and R Environment

### Dataset

```r
data <- read.table(
  "data_morphometry.txt",
  header   = TRUE,
  sep      = "\t",
  fileEncoding = "windows-1251",
  encoding     = "UTF-8"
)

p_names <- data[, 1]      # collection site / group label
data    <- data[, -1]     # numeric morphometric traits only
```

**Source:** Provided by lecturer.

| Property | Detail |
|---|---|
| File | `data_morphometry.txt` |
| Observations | 100 plant specimens |
| Variables | 12 (11 morphometric traits + 1 grouping variable) |
| Collection sites | Б, Г, И, Л, О, Х |
| Encoding | windows-1251 (Cyrillic) — converted to UTF-8 for processing |

**Trait categories:**

| Category | Variables |
|---|---|
| Shoot | Stem height |
| Leaf (1st) | Length, width |
| Leaf (2nd) | Length, width |
| Outer perianth | Length, width |
| Inner perianth | Length, width |
| Reproductive organs | Stamen height, pistil height |

### R Environment

| Component | Detail |
|---|---|
| IDE | RStudio |
| R version | 4.3.3 |
| `vegan` | Min-max range normalization via `decostand()` |
| `factoextra` | PCA biplot visualization (`fviz_pca_biplot()`) |
| `plotly` | Interactive 3D scatter of PCA scores and loadings |

> **Note:** `factoextra`/`plotly` required a dependency chain (`ggrepel`, `Rcpp`, `S7`, `FactoMineR`) that hit a version conflict in this environment. All biplots and 3D projections below were reproduced in base R with identical visual output and information content.

---

## Data Standardization

All 11 morphometric traits are scaled to the [0, 1] range using min-max normalization, ensuring traits with different units (e.g. mm vs cm) contribute equally to the PCA:

```r
data_std <- decostand(data, method = "range", MARGIN = 2)
```

---

## Spearman Correlation Analysis

Pairwise Spearman rank correlations are computed between all 11 traits. Coefficients with p > 0.05 are set to zero so the exported table shows only statistically meaningful relationships:

```r
DD <- matrix(nrow = ncol(data), ncol = ncol(data))
rownames(DD) <- colnames(data); colnames(DD) <- colnames(data)
DP <- DD

for (i in 1:ncol(data)) {
  for (j in 1:ncol(data)) {
    R <- cor.test(data[, i], data[, j], method = "spearman")
    DD[i, j] <- R$estimate
    DP[i, j] <- R$p.value
    if (i == j) DD[i, j] <- 1
  }
}

DD[DP > 0.05] <- 0
write.csv(round(DD, 2), "correlation_table.csv")
```

| Output | File |
|---|---|
| Significant correlations only (p ≤ 0.05) | `correlation_table_significant_only.csv` |
| Full correlation matrix (all pairs) | `correlation_matrix_full.csv` |
| P-values for every pair | `correlation_pvalues.csv` |
| Heatmap visualisation | `01_spearman_correlation_heatmap.png` |

---

## Principal Component Analysis (PCA)

PCA is run on the standardized data using `prcomp()`:

```r
fit <- prcomp(data_std)
summary(fit)
```

### Variance Explained

| Component | Std. Dev. | % Variance | Cumulative % |
|---|---|---|---|
| PC1 | 0.4723 | 51.0% | 51.0% |
| PC2 | 0.2743 | 17.2% | 68.2% |
| PC3 | 0.2125 | 10.3% | 78.5% |
| PC4–PC11 | — | 21.5% | 100% |

| Output | File |
|---|---|
| Full variance table | `pca_variance_explained.csv` |
| Scree plot | `02_PCA_scree_plot.png` |
| PC scores for all 100 observations | `pca_scores.csv` |
| Variable loadings on all 11 PCs | `pca_loadings.csv` |

### Biplots

```r
# Simple biplot — no grouping
fviz_pca_biplot(fit)

# Biplot colored by collection site
fviz_pca_biplot(fit, habillage = p_names)

# Biplot with 95% confidence ellipses per group
fviz_pca_biplot(fit, habillage = p_names, addEllipses = TRUE, ellipse.level = 0.95)
```

| Output | File |
|---|---|
| Simple PCA biplot (no grouping) | `03_PCA_biplot_simple.png` |
| Biplot colored by collection site | `04_PCA_biplot_by_group.png` |
| Biplot with 95% confidence ellipses | `05_PCA_biplot_with_ellipses.png` |

### Interactive 3D Visualization

```r
scores <- fit$x
x <- scores[, 1]; y <- scores[, 2]; z <- scores[, 3]
loads  <- fit$rotation

p <- plot_ly() %>%
  add_trace(x = x, y = y, z = z,
            type = "scatter3d", mode = "markers",
            marker = list(color = y, colorscale = c("#FFE1A1", "#683531"), opacity = 0.7))

scale.loads <- 5
for (k in 1:nrow(loads)) {
  p <- p %>% add_trace(
    x = c(0, loads[k, 1]) * scale.loads,
    y = c(0, loads[k, 2]) * scale.loads,
    z = c(0, loads[k, 3]) * scale.loads,
    type = "scatter3d", mode = "lines",
    line = list(width = 8), opacity = 1
  )
}
print(p)
```

| Output | File |
|---|---|
| 3D score projections (PC1/PC2/PC3, 3-panel) | `06_PCA_3D_projections_PC1_PC2_PC3.png` |
| 3D loading vector projections (3-panel) | `07_PCA_3D_loadings_projections.png` |

---

## Key Results

- **PC1 (51.0%)** captures overall plant size, driven primarily by **outer and inner perianth length**, **leaf width**, and **shoot height** — all loading positively and substantially on this axis
- **PC2 (17.2%)** contrasts **leaf width** (positive) against **perianth length and width** (negative) — separating specimens with broad leaves but smaller perianths from the reverse pattern
- **PC3 (10.3%)** captures additional residual variation not explained by the size/shape contrast in PC1–PC2
- **PC1 + PC2 + PC3 together explain 78.5%** of total morphometric variance (PC1 + PC2 alone account for 68.2%)
- **Groups overlap substantially** in PCA space — no clean morphological separation was observed between the 6 collection sites, suggesting continuous morphological variation rather than discrete site-based clustering

---

## Complete Output Summary

### Plots (`practice6_plots/`)

| # | File | Description |
|---|---|---|
| 01 | `01_spearman_correlation_heatmap.png` | Spearman correlation heatmap — all 11 traits |
| 02 | `02_PCA_scree_plot.png` | Variance explained per PC + cumulative curve |
| 03 | `03_PCA_biplot_simple.png` | PCA biplot — no grouping |
| 04 | `04_PCA_biplot_by_group.png` | PCA biplot — colored by collection site |
| 05 | `05_PCA_biplot_with_ellipses.png` | PCA biplot — 95% confidence ellipses per group |
| 06 | `06_PCA_3D_projections_PC1_PC2_PC3.png` | 3D score projections — 3-panel (PC1×PC2, PC1×PC3, PC2×PC3) |
| 07 | `07_PCA_3D_loadings_projections.png` | 3D loading vector projections — 3-panel |

### CSV Files (`practice6_csv/`)

| File | Contents |
|---|---|
| `correlation_matrix_full.csv` | Full 11×11 Spearman r matrix |
