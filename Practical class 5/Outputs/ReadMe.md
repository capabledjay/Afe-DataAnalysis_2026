# Practice 5: Multivariate Ordination and Clustering of Ecological Community Data

> NMDS Ordination · UPGMA Clustering · envfit Species Vectors · PERMANOVA

---

## Project Overview

This project applies multivariate ordination and clustering techniques to a community ecology dataset. The work progresses through NMDS ordination and UPGMA hierarchical clustering across three distance metrics, species vector fitting with permutation-based significance testing using `envfit`, and PERMANOVA to test whether identified clusters differ significantly in species composition.

---

## Data and R Environment

### Dataset

```r
data <- read.table("data.txt", header = TRUE, sep = "\t", check.names = FALSE)
rownames(data) <- data[, 1]
data <- data[, -1]
```

| Property | Detail |
|---|---|
| File | `data.txt` |
| Sites (rows) | 6 |
| Species (columns) | 34 |
| Value type | Semi-quantitative abundance counts |

**Sites:** Ledyanaya Gora · Karaul village · Ladyginskie Yary · Sopochnaya Karga · Sibiryakov Island Srednee Lake · Chernyi Bay

**Species included:** *Chrysosphaerella brevispina*, *C. coronacircumspina*, *C. longispina*, *Chrysosphaerella* sp., *Paraphysomonas bandaiensis*, *P. gladiata*, *P. vestita*, *Paraphysomonas* sp., *Mallomonas akrokomos*, *M. alata* f. *alata*, *M. alata* f. *hualvensis*, *M. caudata*, *M. crassisquama*, *M. crassisquama* var. *papillosa*, *M. heterospina*, *M. insignis*, *M. striata* var. *striata*, *M. tonsurata*, *Chrysodidymus synuroideus*, *Spiniferomonas abei*, *S. bilacunosa*, *S. bourrellyi*, *S. cornuta*, *S. crusigera*, *S. serrata*, *S. triangularis*, *S. trioralis* f. *trioralis*, *S. trioralis* f. *cuspidata*, *Spiniferomonas* sp., *Synura echinulata* f. *leptorrhabda*, *S. petersenii* f. *petersenii*, *S. petersenii* f. *kufferathii*, *S. punctulosa*, *S. spinosa*

### R Environment

| Component | Detail |
|---|---|
| IDE | RStudio |
| R version | 4.4.0 |
| Key package | `vegan` |

> **Note:** `vegan` is not available via `install.packages()` on some systems due to CRAN restrictions. Install from GitHub source if needed:
> ```r
> curl -L -o vegan.tar.gz "https://github.com/vegandevs/vegan/archive/refs/tags/v2.6-8.tar.gz"
> # then: install.packages("vegan.tar.gz", repos = NULL, type = "source")
> ```

---

## Task 1 — Comparison of Distance Metrics (NMDS and UPGMA)

NMDS ordination and UPGMA hierarchical clustering were performed using three distance metrics. For each metric the same pipeline was applied:

```r
# NMDS ordination
ord <- metaMDS(data, distance = "<metric>")
plot(ord, type = "n")
points(ord, what = "sites", pch = 21, cex = 2.5, lwd = 2.5, col = "<colour>")
text(ord, what = "sites", cex = 0.7, col = "<colour>", pos = 3)

# UPGMA clustering
d   <- vegdist(data, method = "<metric>")
fit <- hclust(d, method = "average")
plot(fit, hang = -1)   # labels aligned at baseline
plot(fit)              # default R style
```

> **Bug fixed:** The original script used `display = "site"` (wrong argument name, wrong value) in all three `text()` calls. The correct argument for `text.ordiplot()` is `what = "sites"`.

### Distance Metrics

| Metric | Function call | Characteristics |
|---|---|---|
| Euclidean | `vegdist(data, method = "euclidean")` | Quantitative distances; sensitive to double zeros |
| Bray–Curtis | `vegdist(data, method = "bray")` | Standard ecological dissimilarity; ignores double zeros |
| Jaccard | `vegdist(data, method = "jaccard")` | Presence/absence-weighted dissimilarity |

**Finding:** All three metrics produced identical clustering topology — a consistent 2-group split was recovered across all methods, indicating a robust community structure signal in the data.

### Output Files (`practice5_plots/`)

| # | File | Description |
|---|---|---|
| 01 | `01_NMDS_euclidean_ordination.png` | NMDS ordination — Euclidean distance |
| 02 | `02_cluster_euclidean_hang.png` | UPGMA dendrogram — Euclidean, labels aligned (`hang = -1`) |
| 03 | `03_cluster_euclidean_default.png` | UPGMA dendrogram — Euclidean, default style |
| 04 | `04_NMDS_bray_curtis_ordination.png` | NMDS ordination — Bray–Curtis distance |
| 05 | `05_cluster_bray_curtis_hang.png` | UPGMA dendrogram — Bray–Curtis, labels aligned |
| 06 | `06_cluster_bray_curtis_default.png` | UPGMA dendrogram — Bray–Curtis, default style |
| 07 | `07_NMDS_jaccard_ordination.png` | NMDS ordination — Jaccard distance |
| 08 | `08_cluster_jaccard_hang.png` | UPGMA dendrogram — Jaccard, labels aligned |
| 09 | `09_cluster_jaccard_default.png` | UPGMA dendrogram — Jaccard, default style |

---

## Task 2 — Detailed Bray–Curtis Analysis

### 2.1 — NMDS (trymax = 50)

```r
ord_bray <- metaMDS(data, distance = "bray", trymax = 50)
```

NMDS stress approached zero due to the small sample size (n = 6 sites). Results are interpreted with caution — near-zero stress with few points indicates the ordination may be fitting the data exactly rather than generalising.

### 2.2 — Significant Species Vectors (envfit)

```r
set.seed(123)
fit_sp <- envfit(ord_bray, data, permutations = 999)

# NOTE: original script used fit_sp$vectors$p.val — correct field is $vectors$pvals
sig_sp <- fit_sp$vectors$pvals <= 0.05
print(fit_sp$vectors$arrows[sig_sp, , drop = FALSE])
print(fit_sp$vectors$pvals[sig_sp])
```

> **Bug fixed:** The original script referenced `fit_sp$vectors$p.val`. The correct field name in `envfit` output is `fit_sp$vectors$pvals`.

**Finding:** Only 1 significant species vector was identified — ***S. cornuta*** (p = 0.033). The low count is expected given the small dataset (6 sites, 720 possible permutations — full enumeration was used automatically).

```
Significant species (p <= 0.05):
              NMDS1      NMDS2
S. cornuta  -0.6641  -0.7476

p-value: 0.0333
```

### 2.3 — UPGMA Clustering (Bray–Curtis)

```r
d_bray  <- vegdist(data, method = "bray")
hc_bray <- hclust(d_bray, method = "average")
```

### 2.4 — Cluster Extraction (cutree, k = 2)

```r
k_clusters <- 2
clusters   <- cutree(hc_bray, k = k_clusters)
clusters   <- factor(clusters, labels = c("Cluster 1", "Cluster 2"))
```

| Cluster | Sites |
|---|---|
| Cluster 1 (blue) | Ledyanaya Gora · Ladyginskie Yary · Sopochnaya Karga · Chernyi Bay |
| Cluster 2 (red) | Karaul village · Sibiryakov Island Srednee Lake |

### 2.5 — Full NMDS Plot with Ellipses and Vectors

```r
cols <- c("blue", "red")[as.numeric(clusters)]

plot(ord_bray, type = "n",
     main = paste0("NMDS (Bray-Curtis), stress = ", round(ord_bray$stress, 4)),
     sub  = "Clusters assigned by UPGMA (k = 2)")

points(ord_bray, col = cols, pch = 16, cex = 2.5)

# 95% SD confidence ellipses per cluster
ordiellipse(ord_bray, groups = clusters,
            col = c("blue", "red"), kind = "sd", conf = 0.95, lwd = 2)

# Non-overlapping site labels
orditorp(ord_bray, display = "sites", cex = 0.8, col = "black", air = 0.5)

# Significant species vectors (p <= 0.05)
plot(fit_sp, p.max = 0.05, col = "darkgreen", cex = 0.9, add = TRUE)

legend("topright", legend = levels(clusters),
       col = c("blue", "red"), pch = 16, title = "Cluster", bty = "n")
```

> **Note:** Cluster 2 contains only 2 sites, so its ellipse is drawn as a line rather than a closed ellipse (rank-deficient covariance matrix). This is expected for groups with n < 3.

### Task 2 Output Files

| # | File | Description |
|---|---|---|
| 10 | `10_NMDS_bray_clusters_ellipses_vectors.png` | Full NMDS — cluster colours, 95% SD ellipses, site labels, *S. cornuta* vector |
| 11 | `11_NMDS_bray_envfit_significant_vectors.png` | NMDS + significant species vectors only (clean view) |
| 12 | `12_cluster_bray_coloured_by_group.png` | UPGMA dendrogram with cluster rectangles (blue = Cluster 1, red = Cluster 2) |
| 14 | `14_envfit_species_pvalues_ranked.png` | All 34 species p-values ranked; green = significant (p ≤ 0.05) |

---

## Task 3 — PERMANOVA

PERMANOVA tests whether the two UPGMA clusters differ significantly in species composition using Bray–Curtis dissimilarities.

```r
set.seed(456)
adonis_result <- adonis2(data ~ clusters, method = "bray", permutations = 999)
print(adonis_result)

R2    <- adonis_result$R2[1]
p_val <- adonis_result$`Pr(>F)`[1]
```

> **Note:** With only 6 sites, only 720 unique permutations are possible. `adonis2` switches to full enumeration automatically, reducing the effective permutation count to 719.

### Results

| Source | Df | SumOfSqs | R² | F | Pr(>F) |
|---|---|---|---|---|---|
| Model (clusters) | 1 | 0.49182 | 0.4088 | 2.7661 | 0.0667 |
| Residual | 4 | 0.71121 | 0.5912 | | |
| Total | 5 | 1.20303 | 1.0000 | | |

```
=======================================
PERMANOVA results (Bray-Curtis):
R² = 0.4088
p-value = 0.06667
Conclusion: differences between clusters are not significant (p >= 0.05).
=======================================

Cluster composition:
Cluster 1: Ledyanaya Gora  Ladyginskie Yary  Sopochnaya Karga  Chernyi Bay
Cluster 2: Karaul village  Sibiryakov Island, Srednee Lake
```

**Interpretation:**
- **R² = 0.4088** — cluster membership explains ~41% of the total variation in species composition
- **p = 0.0667** — the observed separation is not significant at the α = 0.05 threshold
- The lack of significance is most likely a **statistical power issue**: with only 6 sites and 720 possible permutations, the test has very limited ability to detect true differences even when effect sizes are moderate

### Task 3 Output Files

| # | File | Description |
|---|---|---|
| 13 | `13_PERMANOVA_R2_summary.png` | R² bar chart — variance explained by clusters vs residual |

---

## Complete Output Summary

All plots are saved to `practice5_plots/`:

| # | File | Task | Description |
|---|---|---|---|
| 01 | `01_NMDS_euclidean_ordination.png` | 1 | NMDS — Euclidean |
| 02 | `02_cluster_euclidean_hang.png` | 1 | Dendrogram — Euclidean, aligned |
| 03 | `03_cluster_euclidean_default.png` | 1 | Dendrogram — Euclidean, default |
| 04 | `04_NMDS_bray_curtis_ordination.png` | 1 | NMDS — Bray–Curtis |
| 05 | `05_cluster_bray_curtis_hang.png` | 1 | Dendrogram — Bray–Curtis, aligned |
| 06 | `06_cluster_bray_curtis_default.png` | 1 | Dendrogram — Bray–Curtis, default |
| 07 | `07_NMDS_jaccard_ordination.png` | 1 | NMDS — Jaccard |
| 08 | `08_cluster_jaccard_hang.png` | 1 | Dendrogram — Jaccard, aligned |
| 09 | `09_cluster_jaccard_default.png` | 1 | Dendrogram — Jaccard, default |
| 10 | `10_NMDS_bray_clusters_ellipses_vectors.png` | 2.5 | Full NMDS — clusters + ellipses + *S. cornuta* vector |
| 11 | `11_NMDS_bray_envfit_significant_vectors.png` | 2.2 | NMDS — significant envfit vectors only |
| 12 | `12_cluster_bray_coloured_by_group.png` | 2.3 | UPGMA dendrogram with cluster rectangles |
| 13 | `13_PERMANOVA_R2_summary.png` | 3 | PERMANOVA R² bar chart |
| 14 | `14_envfit_species_pvalues_ranked.png` | 2.2 | All 34 species p-values ranked |
