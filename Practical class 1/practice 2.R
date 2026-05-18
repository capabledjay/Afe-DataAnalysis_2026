#--------------------start-------------------------------
# Get current working directory
getwd()
#----------------read dataset--------------------------
example_df<-read.csv("C:/Users/DELL/Desktop/Year_1/data analysis/Practical class 2/distribution.csv", header = TRUE,dec = ',', sep = ";")
factor_df <- read.csv("C:/Users/DELL/Desktop/Year_1/data analysis/Practical class 2/factor_data.csv")
imputed_df <- read.csv("C:/Users/DELL/Desktop/Year_1/data analysis/Practical class 2/imputed_data.csv")
# Display structure with variable types
str(example_df)
str(factor_df)
str(imputed_df)
#---------------merge two files-------------------------
data_for_analysis <- merge(
  factor_df, 
  imputed_df, 
  by = "record_id",        # column for merge
  all = FALSE       # FALSE = INNER JOIN (only coincidences), TRUE = FULL JOIN
)
str(data_for_analysis)

# save data_for_analysis in CSV
write.csv(data_for_analysis, "data_for_analysis.csv", row.names = FALSE)  

#------------------Probability Distributions----------------------- 
install.packages("MASS", dependencies=T)
library(MASS)
#----------------example-----------------------------------------
summary (example_df)
example_df$value <- as.numeric(example_df$value)
summary(example_df)
#building histograms for example
# normal distribution
val<-example_df[example_df$distribution=="norm",]$value

mean(val)

sd(val)

hist(val)

fit<-fitdistr(val, densfun="normal")

fit
#lognormal distribution
val<-example_df[example_df$distribution=="lognorm",]$value

mean(val)

sd(val)

hist(val)

fit<-fitdistr(val, densfun="lognormal")

fit

unname(fit$estimate[1])

unname(fit$estimate[2])

m_log<-exp(unname(fit$estimate[1]))*sqrt(exp(unname(fit$estimate[2])^2))
m_log
sd_log<-sqrt(exp(2*unname(fit$estimate[1]))*(exp(unname(fit$estimate[2])^2)-1)*sqrt(exp(unname(fit$estimate[2])^2)))
sd_log
#exponential distribution

val<-example_df[example_df$distribution=="exp",]$value

mean(val)

sd(val)

hist(val)

fit<-fitdistr(val, densfun="exponential")

fit

unname(fit$estimate[1])

m_exp<-1/unname(fit$estimate[1])
m_exp
#Poisson distribution
val<-example_df[example_df$distribution=="pois",]$value

mean(val)

sd(val)

hist(val)

fit<-fitdistr(val, densfun="Poisson")

fit

unname(fit$estimate[1])

sd_pois<-sqrt(unname(fit$estimate[1]))
sd_pois
#Selecting a Distribution Model

val<-example_df[example_df$distribution=="lognorm",]$value

fit_1<-fitdistr(val, densfun="normal")
fit_2<-fitdistr(val, densfun="lognormal")
fit_3<-fitdistr(val, densfun="exponential")

#Bayesian Information Criterion calculation
BIC(fit_3)

#calculation of the Bayesian information criterion for all models
BIC_value<-c(BIC(fit_1), BIC(fit_2), BIC(fit_3))

#forming a vector with the name of the models
distribution<-c("normal", "lognormal", "exponential")

#combining the results into a final table
rez<-data.frame(BIC_value=BIC_value, distribution=distribution)

#sort table in ascending order of Bayesian Information Criterion value
rez<-rez[order(rez$BIC_value, decreasing=F),]

rez


#calculation of absolute values of the confidence interval for the mean of a lognormal distribution
error_min<-unname(fit_2$estimate[1])-unname(fit_2$sd[1])
error_max<-unname(fit_2$estimate[1])+unname(fit_2$sd[1])

error_min
error_max

m<-exp(unname(fit_2$estimate[1]))*sqrt(exp(unname(fit_2$estimate[2])^2))
value_error_min<-exp(error_min)*sqrt(exp(unname(fit_2$estimate[2])^2))
value_error_max<-exp(error_max)*sqrt(exp(unname(fit_2$estimate[2])^2))

value_error_min
m
value_error_max

#--------------data for analysis--------------------------
#building histograms
value_d1<-data_for_analysis$lipids1
hist(value_d1)
value_d2<-data_for_analysis$lipids2
hist(value_d2)
value_d3<-data_for_analysis$lipids3
hist(value_d3)
value_d4<-data_for_analysis$lipids4
hist(value_d4)


# d1 distribution estimate


fit_d1_1<-fitdistr(value_d1,densfun="normal")
fit_d1_2<-fitdistr(value_d1,densfun="lognormal")
fit_d1_3<-fitdistr(value_d1,densfun="exponential")

#calculation of the Bayesian information criterion (BIC) and finding of BIC minimum for d1

BIC_value_d1 <- c(BIC(fit_d1_1),BIC(fit_d1_2),BIC(fit_d1_3))
distribution <-c("normal","lognormal","exponential")
result_d1<-data.frame(BIC_value_d1=BIC_value_d1, distribution=distribution)
result_d1
min(result_d1$BIC_value_d1)
distribution_d1<-result_d1[result_d1$BIC_value_d1==min(result_d1$BIC_value_d1),]$distribution
distribution_d1
# Finding parameters for d1
fit_d1_1$estimate[1:2]

# d2 distribution estimate


fit_d2_1<-fitdistr(value_d2,densfun="normal")
fit_d2_2<-fitdistr(value_d2,densfun="lognormal")
fit_d2_3<-fitdistr(value_d2,densfun="exponential")

#calculation of the Bayesian information criterion (BIC) and finding of BIC minimum for d2

BIC_value_d2 <- c(BIC(fit_d2_1),BIC(fit_d2_2),BIC(fit_d2_3))
distribution <-c("normal","lognormal","exponential")
result_d2<-data.frame(BIC_value_d2=BIC_value_d2, distribution=distribution)
result_d2
min(result_d2$BIC_value_d2)
distribution_d2<-result_d2[result_d2$BIC_value_d2==min(result_d2$BIC_value_d2),]$distribution
distribution_d2
# Finding parameters for d2
fit_d2_1$estimate[1:2]

#================================================================
# TASK 1: Distribution Analysis by Outcome Group
# Goal: Identify best-fitting distribution for each continuous
# variable by outcome group, then compile descriptive statistics
#================================================================

#----------------------------------------------------------------
# SECTION 1: Data Preparation
# Split data by outcome: 0 = no tumour, 1 = tumour
#----------------------------------------------------------------

group0 <- data_for_analysis[data_for_analysis$outcome == 0, ]
group1 <- data_for_analysis[data_for_analysis$outcome == 1, ]

# Expected: Total=1148, Group0=988, Group1=161
# Note: 1 observation has NA outcome, excluded from both groups
cat("Group sizes - Total:", nrow(data_for_analysis), 
    "| G0:", nrow(group0), 
    "| G1:", nrow(group1), "\n")

# Select continuous variables by column index to avoid encoding issues
# Column 20 (lipids5) excluded per standard task instructions
continuous_vars <- colnames(data_for_analysis)[c(
  7,   # hormone1
  8,   # hormone2
  9,   # hormone3
  10,  # hormone4
  11,  # hormone5
  12,  # hormone6
  13,  # hormone7
  14,  # hormone8
  15,  # hormone10_generated
  16,  # lipids1
  17,  # lipids2
  18,  # lipids3
  19,  # lipids4
  21,  # carb_metabolism (Cyrillic "c" - accessed by index)
  22,  # lipid_pero1
  23,  # lipid_pero2
  24,  # lipid_pero3
  25,  # lipid_pero4
  26,  # lipid_pero5
  27,  # antioxidant1
  28,  # antioxidant2
  29,  # antioxidant3
  30,  # antioxidant4
  31   # antioxidant5
)]
# Expected: 24 variables
cat("Variables selected:", length(continuous_vars), "\n")


#----------------------------------------------------------------
# SECTION 2: Distribution Fitting Function
# Fits normal, lognormal, exponential to a data vector
# Selects best fit by lowest BIC and returns its parameters
# Note: lognormal/exponential require strictly positive values
#----------------------------------------------------------------

fit_best_distribution <- function(x, var_name, group_name) {
  
  x <- x[!is.na(x)]
  
  if (length(x) < 5) {
    return(data.frame(
      variable = var_name, group = group_name, n = length(x),
      best_dist = NA, param1_name = NA, param1_value = NA,
      param2_name = NA, param2_value = NA,
      BIC_normal = NA, BIC_lognormal = NA, BIC_exponential = NA
    ))
  }
  
  all_positive <- all(x > 0)
  fit_norm  <- fitdistr(x, densfun = "normal")
  bic_norm  <- BIC(fit_norm)
  
  if (all_positive) {
    fit_lognorm <- fitdistr(x, densfun = "lognormal")
    fit_exp     <- fitdistr(x, densfun = "exponential")
    bic_lognorm <- BIC(fit_lognorm)
    bic_exp     <- BIC(fit_exp)
  } else {
    # Zero/negative values present - only normal applicable
    bic_lognorm <- NA
    bic_exp     <- NA
  }
  
  bic_values <- c(normal = bic_norm, lognormal = bic_lognorm, exponential = bic_exp)
  best_dist  <- names(which.min(bic_values))
  
  if (best_dist == "normal") {
    param1_name  <- "mean"
    param1_value <- round(unname(fit_norm$estimate["mean"]), 4)
    param2_name  <- "sd"
    param2_value <- round(unname(fit_norm$estimate["sd"]), 4)
  } else if (best_dist == "lognormal") {
    param1_name  <- "meanlog"
    param1_value <- round(unname(fit_lognorm$estimate["meanlog"]), 4)
    param2_name  <- "sdlog"
    param2_value <- round(unname(fit_lognorm$estimate["sdlog"]), 4)
  } else if (best_dist == "exponential") {
    param1_name  <- "rate"
    param1_value <- round(unname(fit_exp$estimate["rate"]), 4)
    param2_name  <- "mean(1/rate)"
    param2_value <- round(1 / unname(fit_exp$estimate["rate"]), 4)
  }
  
  return(data.frame(
    variable = var_name, group = group_name, n = length(x),
    best_dist = best_dist,
    param1_name = param1_name, param1_value = param1_value,
    param2_name = param2_name, param2_value = param2_value,
    BIC_normal = round(bic_norm, 4),
    BIC_lognormal = round(bic_lognorm, 4),
    BIC_exponential = round(bic_exp, 4)
  ))
}


#----------------------------------------------------------------
# SECTION 3: Run Distribution Fitting for All Variables
# Loops through all 24 variables x 2 groups = 48 combinations
# Expected output: results_table with 48 rows
#----------------------------------------------------------------

results_list <- list()
counter <- 1

for (var in continuous_vars) {
  results_list[[counter]]   <- fit_best_distribution(group0[[var]], var, "0 - no tumour")
  results_list[[counter+1]] <- fit_best_distribution(group1[[var]], var, "1 - tumour")
  counter <- counter + 2
}

results_table <- do.call(rbind, results_list)
rownames(results_table) <- NULL
# Expected: 48 rows
cat("Distribution fitting complete. Rows:", nrow(results_table), "\n")

# REQUIREMENT 1 RESULT: Best fit distributions by outcome group
# Expected: most variables lognormal, some normal where zeros present
print(results_table[, c("variable", "group", "best_dist",
                        "BIC_normal", "BIC_lognormal", "BIC_exponential")])

# Distribution frequency per group
# Expected G0: 18 lognormal, 6 normal | G1: 20 lognormal, 4 normal
cat("\nG0 distribution frequency:\n")
print(table(results_table$best_dist[results_table$group == "0 - no tumour"]))
cat("G1 distribution frequency:\n")
print(table(results_table$best_dist[results_table$group == "1 - tumour"]))

# Variables where distribution differs between groups
# Expected: 6 variables differ (hormone3, hormone7, antioxidant3,
# antioxidant4, lipid_pero5, carb_metabolism)
g0 <- results_table[results_table$group == "0 - no tumour", c("variable", "best_dist")]
g1 <- results_table[results_table$group == "1 - tumour",    c("variable", "best_dist")]
diff_vars <- merge(g0, g1, by = "variable", suffixes = c("_g0", "_g1"))
diff_vars <- diff_vars[diff_vars$best_dist_g0 != diff_vars$best_dist_g1, ]
cat("Variables with different distributions between groups:\n")
print(diff_vars)


#----------------------------------------------------------------
# SECTION 4: Descriptive Statistics Table
# Combines basic stats with distribution parameters from Section 3
# REQUIREMENT 2 RESULT: saved as descriptive_statistics_by_group.csv
#----------------------------------------------------------------

get_descriptive_stats <- function(x, var_name, group_name) {
  x <- x[!is.na(x)]
  data.frame(
    variable = var_name, group = group_name,
    n        = length(x),
    mean     = round(mean(x), 4),
    sd       = round(sd(x), 4),
    median   = round(median(x), 4),
    min      = round(min(x), 4),
    max      = round(max(x), 4),
    q25      = round(quantile(x, 0.25), 4),
    q75      = round(quantile(x, 0.75), 4)
  )
}

desc_list    <- list()
desc_counter <- 1

for (var in continuous_vars) {
  desc_list[[desc_counter]]   <- get_descriptive_stats(group0[[var]], var, "0 - no tumour")
  desc_list[[desc_counter+1]] <- get_descriptive_stats(group1[[var]], var, "1 - tumour")
  desc_counter <- desc_counter + 2
}

desc_table <- do.call(rbind, desc_list)
rownames(desc_table) <- NULL

# Merge descriptive stats with distribution results
final_table <- merge(
  desc_table,
  results_table[, c("variable", "group", "best_dist",
                    "param1_name", "param1_value",
                    "param2_name", "param2_value",
                    "BIC_normal", "BIC_lognormal", "BIC_exponential")],
  by = c("variable", "group")
)

final_table <- final_table[order(final_table$variable, final_table$group), ]
rownames(final_table) <- NULL

# Save final table - main deliverable for GitHub submission
# Expected: 48 rows, 17 columns
write.csv(final_table, "descriptive_statistics_by_group.csv", row.names = FALSE)
cat("Saved: descriptive_statistics_by_group.csv |",
    nrow(final_table), "rows |", ncol(final_table), "cols\n")

# Print final table to console
print(final_table[, c("variable", "group", "n",
                      "mean", "sd", "median",
                      "min", "max", "q25", "q75",
                      "best_dist",
                      "param1_name", "param1_value",
                      "param2_name", "param2_value")])



#================================================================
# EXTRA POINTS TASK: Find and Fix Error in lipids5
#================================================================

#----------------------------------------------------------------
# SECTION 5: Investigate lipids5
# Explore the variable to identify the deliberate error
# Expected: summary stats, NA count, histogram
#----------------------------------------------------------------

lipids5_vals <- data_for_analysis$lipids5

# Summary and NA check
# Expected: 276 NAs, values range 0.09-1.24
print(summary(lipids5_vals))
cat("Total NAs:", sum(is.na(lipids5_vals)), 
    "| Non-NA:", sum(!is.na(lipids5_vals)), "\n")

# Compare with lipids1-4 to check if missingness is unique to lipids5
# Expected: lipids1-4 have no NAs, only lipids5 has 276 NAs
print(summary(data_for_analysis[, c("lipids1","lipids2",
                                    "lipids3","lipids4",
                                    "lipids5")]))

# Histogram to visually inspect distribution shape
# Expected: right-skewed distribution, no extreme outliers
hist(lipids5_vals,
     main = "Histogram of lipids5 - checking for errors",
     xlab = "lipids5 values",
     col  = "lightblue")

#----------------------------------------------------------------
# SECTION 5b: Investigate the missing data pattern in lipids5
# Expected: NAs scattered across dataset, disproportionate
# between groups - indicating systematic missingness (MNAR)
#----------------------------------------------------------------

# NA count by outcome group
# Expected: G0=248/988, G1=29/161 - disproportionate missingness
cat("lipids5 NAs by group - G0:", sum(is.na(group0$lipids5)),
    "of", nrow(group0),
    "| G1:", sum(is.na(group1$lipids5)),
    "of", nrow(group1), "\n")

# Check if NAs are isolated to lipids5 or affect other lipids too
# Expected: 0 NAs in lipids1-4 for the same missing records
# confirming lipids5 was specifically missed during imputation
missing_lipids5_rows <- data_for_analysis[is.na(data_for_analysis$lipids5), ]
cat("NAs in other lipids for same records - lipids1:",
    sum(is.na(missing_lipids5_rows$lipids1)),
    "| lipids2:", sum(is.na(missing_lipids5_rows$lipids2)),
    "| lipids3:", sum(is.na(missing_lipids5_rows$lipids3)),
    "| lipids4:", sum(is.na(missing_lipids5_rows$lipids4)), "\n")

# Check spread of missing record IDs
# Expected: NAs scattered from record 13 to 2694 (full range)
# ruling out a batch entry error
na_records <- data_for_analysis$record_id[is.na(data_for_analysis$lipids5)]
cat("Missing record_id range:", min(na_records), "to", max(na_records),
    "| Full dataset range:", min(data_for_analysis$record_id),
    "to", max(data_for_analysis$record_id), "\n")

# CONCLUSION: lipids5 was omitted from the imputation pipeline in
# Practice 1. All other lipids were fully imputed while lipids5
# was left with 276 NAs (24% missing) scattered across the full
# dataset. Missingness is disproportionate between groups (25% G0,
# 18% G1) indicating MNAR - systematic rather than random missingness.

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
# SECTION 5c: Trace the source of missing lipids5
# Verify whether the NAs originated in the source data
# Expected: 276 NAs already present in imputed_df, confirming
# lipids5 was never imputed in Practice 1
#----------------------------------------------------------------

# Expected: same 276 NAs as data_for_analysis - error is in source
cat("lipids5 NAs in imputed_df:", sum(is.na(imputed_df$lipids5)),
    "| in data_for_analysis:", sum(is.na(data_for_analysis$lipids5)), "\n")

# Expected: lipids5 not in factor_df - only exists in imputed_df
cat("lipids5 in factor_df:", "lipids5" %in% colnames(factor_df),
    "| in imputed_df:", "lipids5" %in% colnames(imputed_df), "\n")

# CONCLUSION: 276 NAs were already present in imputed_df confirming
# lipids5 was skipped during the Practice 1 imputation pipeline.
# The error originated in the source imputed dataset, not the merge.

#----------------------------------------------------------------
# SECTION 5d: Fix lipids5 - Random Sampling Imputation
# Median imputation was rejected as it created an artificial spike
# at 0.41 destroying the natural distribution shape.
# Random sampling from observed values within each group preserves
# the original right-skewed distribution shape.
# set.seed(123) used for reproducibility.
#----------------------------------------------------------------

data_for_analysis_fixed <- data_for_analysis
set.seed(123)

# Identify NA positions by group explicitly
# Using which() avoids issues with NA outcome comparisons
na_pos_g0 <- which(is.na(data_for_analysis_fixed$lipids5) &
                     !is.na(data_for_analysis_fixed$outcome) &
                     data_for_analysis_fixed$outcome == 0)

na_pos_g1 <- which(is.na(data_for_analysis_fixed$lipids5) &
                     !is.na(data_for_analysis_fixed$outcome) &
                     data_for_analysis_fixed$outcome == 1)

na_pos_na_outcome <- which(is.na(data_for_analysis_fixed$lipids5) &
                             is.na(data_for_analysis_fixed$outcome))

# Expected: G0=247, G1=28, NA outcome=1
cat("NAs to impute - G0:", length(na_pos_g0),
    "| G1:", length(na_pos_g1),
    "| NA outcome:", length(na_pos_na_outcome), "\n")

# Get observed values per group for sampling
observed_g0  <- data_for_analysis_fixed$lipids5[
  !is.na(data_for_analysis_fixed$lipids5) &
    !is.na(data_for_analysis_fixed$outcome) &
    data_for_analysis_fixed$outcome == 0]

observed_g1  <- data_for_analysis_fixed$lipids5[
  !is.na(data_for_analysis_fixed$lipids5) &
    !is.na(data_for_analysis_fixed$outcome) &
    data_for_analysis_fixed$outcome == 1]

observed_all <- data_for_analysis_fixed$lipids5[
  !is.na(data_for_analysis_fixed$lipids5)]

# Expected: G0=740, G1=132 observed values available
cat("Observed values - G0:", length(observed_g0),
    "| G1:", length(observed_g1), "\n")

# Impute by random sampling within group
data_for_analysis_fixed$lipids5[na_pos_g0] <- sample(
  observed_g0, size = length(na_pos_g0), replace = TRUE)

data_for_analysis_fixed$lipids5[na_pos_g1] <- sample(
  observed_g1, size = length(na_pos_g1), replace = TRUE)

if (length(na_pos_na_outcome) > 0) {
  data_for_analysis_fixed$lipids5[na_pos_na_outcome] <- sample(
    observed_all, size = length(na_pos_na_outcome), replace = TRUE)
}

# Verify fix - Expected: 0 NAs remaining
# Summary stats should be virtually unchanged from original
cat("NAs remaining:", sum(is.na(data_for_analysis_fixed$lipids5)), "\n")
cat("Before fix:\n")
print(summary(data_for_analysis$lipids5))
cat("After fix:\n")
print(summary(data_for_analysis_fixed$lipids5))

# Compare histograms - After Fix shape should match Before Fix
# Expected: same right-skewed shape preserved, no artificial spikes
par(mfrow = c(1, 2))
hist(data_for_analysis$lipids5,
     main = "lipids5 - Before Fix",
     xlab = "lipids5", col = "lightblue", breaks = 20)
hist(data_for_analysis_fixed$lipids5,
     main = "lipids5 - After Fix (Random Sample)",
     xlab = "lipids5", col = "lightgreen", breaks = 20)
par(mfrow = c(1, 1))

# Save fixed dataset
# Expected: data_for_analysis_fixed.csv with 0 NAs in lipids5
write.csv(data_for_analysis_fixed,
          "data_for_analysis_fixed.csv",
          row.names = FALSE)
cat("Saved: data_for_analysis_fixed.csv\n")


#----------------------------------------------------------------
# SECTION 6: Re-run Distribution Analysis with lipids5 included
# Using data_for_analysis_fixed with imputed lipids5
# Reuses fit_best_distribution() from Section 2
#----------------------------------------------------------------

# Update groups to use fixed dataset
group0_fixed <- data_for_analysis_fixed[data_for_analysis_fixed$outcome == 0, ]
group1_fixed <- data_for_analysis_fixed[data_for_analysis_fixed$outcome == 1, ]

# Include lipids5 (col 20) - now fixed
# Expected: 25 variables
continuous_vars_fixed <- colnames(data_for_analysis_fixed)[c(
  7,   # hormone1
  8,   # hormone2
  9,   # hormone3
  10,  # hormone4
  11,  # hormone5
  12,  # hormone6
  13,  # hormone7
  14,  # hormone8
  15,  # hormone10_generated
  16,  # lipids1
  17,  # lipids2
  18,  # lipids3
  19,  # lipids4
  20,  # lipids5 - included (fixed)
  21,  # carb_metabolism
  22,  # lipid_pero1
  23,  # lipid_pero2
  24,  # lipid_pero3
  25,  # lipid_pero4
  26,  # lipid_pero5
  27,  # antioxidant1
  28,  # antioxidant2
  29,  # antioxidant3
  30,  # antioxidant4
  31   # antioxidant5
)]
cat("Variables selected:", length(continuous_vars_fixed), "\n")

# Run distribution fitting loop
# Expected: 50 rows (25 variables x 2 groups)
results_list_fixed <- list()
counter_fixed <- 1

for (var in continuous_vars_fixed) {
  results_list_fixed[[counter_fixed]]   <- fit_best_distribution(
    group0_fixed[[var]], var, "0 - no tumour")
  results_list_fixed[[counter_fixed+1]] <- fit_best_distribution(
    group1_fixed[[var]], var, "1 - tumour")
  counter_fixed <- counter_fixed + 2
}

results_table_fixed <- do.call(rbind, results_list_fixed)
rownames(results_table_fixed) <- NULL
# Expected: 50 rows
cat("Rows:", nrow(results_table_fixed), "\n")

# REQUIREMENT 1 RESULT (extra points version)
# Expected: most variables lognormal, lipids5 lognormal both groups
print(results_table_fixed[, c("variable", "group", "best_dist",
                              "BIC_normal", "BIC_lognormal",
                              "BIC_exponential")])

# Distribution frequency per group
# Expected G0: 19 lognormal, 6 normal | G1: 21 lognormal, 4 normal
cat("\nG0 distribution frequency:\n")
print(table(results_table_fixed$best_dist[
  results_table_fixed$group == "0 - no tumour"]))
cat("G1 distribution frequency:\n")
print(table(results_table_fixed$best_dist[
  results_table_fixed$group == "1 - tumour"]))

# Variables where distribution differs between groups
# Expected: same 6 variables as standard task
g0_fixed <- results_table_fixed[
  results_table_fixed$group == "0 - no tumour", c("variable", "best_dist")]
g1_fixed <- results_table_fixed[
  results_table_fixed$group == "1 - tumour", c("variable", "best_dist")]
diff_vars_fixed <- merge(g0_fixed, g1_fixed, by = "variable",
                         suffixes = c("_g0", "_g1"))
diff_vars_fixed <- diff_vars_fixed[
  diff_vars_fixed$best_dist_g0 != diff_vars_fixed$best_dist_g1, ]
cat("Variables with different distributions between groups:\n")
print(diff_vars_fixed)

#----------------------------------------------------------------
# SECTION 7: Final Descriptive Statistics Table (extra points)
# Reuses get_descriptive_stats() from Section 4
# Expected: 50 rows, 18 columns
# Saved as: descriptive_statistics_by_groups_fixed.csv
#----------------------------------------------------------------

desc_list_fixed    <- list()
desc_counter_fixed <- 1

for (var in continuous_vars_fixed) {
  desc_list_fixed[[desc_counter_fixed]]   <- get_descriptive_stats(
    group0_fixed[[var]], var, "0 - no tumour")
  desc_list_fixed[[desc_counter_fixed+1]] <- get_descriptive_stats(
    group1_fixed[[var]], var, "1 - tumour")
  desc_counter_fixed <- desc_counter_fixed + 2
}

desc_table_fixed <- do.call(rbind, desc_list_fixed)
rownames(desc_table_fixed) <- NULL

# Merge descriptive stats with distribution results
final_table_fixed <- merge(
  desc_table_fixed,
  results_table_fixed[, c("variable", "group", "best_dist",
                          "param1_name", "param1_value",
                          "param2_name", "param2_value",
                          "BIC_normal", "BIC_lognormal",
                          "BIC_exponential")],
  by = c("variable", "group")
)

final_table_fixed <- final_table_fixed[
  order(final_table_fixed$variable, final_table_fixed$group), ]
rownames(final_table_fixed) <- NULL

# Save - main extra points deliverable for GitHub
write.csv(final_table_fixed,
          "descriptive_statistics_by_groups_fixed.csv",
          row.names = FALSE)
cat("Saved: descriptive_statistics_by_groups_fixed.csv |",
    nrow(final_table_fixed), "rows |", ncol(final_table_fixed), "cols\n")

# Print final table
print(final_table_fixed[, c("variable", "group", "n",
                            "mean", "sd", "median",
                            "min", "max", "q25", "q75",
                            "best_dist",
                            "param1_name", "param1_value",
                            "param2_name", "param2_value")])

# EXTRA POINTS TASK - SUMMARY
# Error:    lipids5 had 276 NAs (24%) while lipids1-4 had none
#           confirming it was skipped in Practice 1 imputation
# Fix:      Random sampling from observed values within each group
#           G0: 247 NAs filled from 740 observed values
#           G1:  28 NAs filled from 132 observed values
#           NA outcome: 1 NA filled from all observed values
# Result:   0 NAs remaining, distribution shape preserved
#           Mean: 0.4333 -> 0.4287 | Q1/Q3 unchanged
#           lipids5 follows lognormal in both groups
# Outputs:  data_for_analysis_fixed.csv
#           descriptive_statistics_by_groups_fixed.csv

