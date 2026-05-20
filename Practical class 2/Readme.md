PROJECT DESCRIPTION

This project applies statistical distribution analysis to a clinical dataset. The main work involved splitting data by outcome group, fitting probability distributions using Maximum Likelihood Estimation, selecting the best fit using the Bayesian Information Criterion, and building a descriptive statistics table by group. An extra points task required finding and fixing a missing data error in the lipids5 variable before repeating the analysis.

DATA DESCRIPTION

The dataset "data_for_analysis" was built by merging two files from Practice 1 — factor_df containing categorical and outcome variables, and imputed_handle_MD_df_final containing the imputed continuous variables. The final merged dataset has 1,148 records and 31 variables covering hormone measurements (hormone1 to hormone8, hormone10_generated), lipid markers (lipids1 to lipids5, lipid_pero1 to lipid_pero5), antioxidant indices (antioxidant1 to antioxidant5), carbohydrate metabolism, and factor variables (outcome, factor_eth, factor_h, factor_pcos, factor_prl). The outcome variable (0 = no tumour, 1 = tumour) was used as the grouping variable throughout.

R Environment

IDE: RStudio (Mac), R version 4.5.3

Key Packages Used

MASS — for fitdistr() to fit probability distributions using MLE and calculate BIC.

PROCEDURES

The dataset was split by outcome into Group 0 (no tumour, 988 observations) and Group 1 (tumour, 161 observations). One observation with a missing outcome value was excluded from both groups. Continuous variables were selected by column index rather than name to avoid an encoding issue with the carb_metabolism column, which had a Cyrillic character in its name. A total of 24 continuous variables were used for the standard task, with lipids5 excluded per task instructions.


A reusable function fit_best_distribution() was written to handle the fitting process. For each variable and group, it removes NAs, fits normal, lognormal, and exponential distributions, compares BIC values, picks the lowest, and returns the parameters of the winning distribution. Variables with zero or negative values were restricted to normal distribution fitting only since lognormal and exponential require strictly positive inputs.

TASK 1 — STANDARD VERSION

The function was applied across all 24 variables and both outcome groups, giving 48 results in total.

In Group 0, 18 variables followed a lognormal distribution and 6 followed normal. In Group 1, 20 were lognormal and 4 were normal. The normal fits occurred in variables containing zeros — hormone3, hormone5, hormone7, and antioxidant4. Six variables showed different distributions between the two groups: hormone3, hormone7, antioxidant3, and antioxidant4 shifted from normal in Group 0 to lognormal in Group 1, while lipid_pero5 and carb_metabolism went the other way. This likely reflects biological differences associated with tumour presence.

The final descriptive statistics table combines n, mean, sd, median, min, max, Q25, Q75, the best-fit distribution, and its parameters for each variable and group. It contains 48 rows and 17 columns.

Output file: descriptive_statistics_by_group.csv

EXTRA POINTS — FINDING AND FIXING THE ERROR IN lipids5

Investigation started with a histogram and summary statistics. The distribution shape looked reasonable — right-skewed, values between 0.09 and 1.24, no extreme outliers. The problem was 276 missing values, 24% of the dataset. What made this clearly an error was that lipids1 through lipids4 had zero NAs in those exact same records, meaning lipids5 was the only lipid variable skipped during the Practice 1 imputation. The missingness was also uneven between groups — 25% in Group 0 versus 18% in Group 1 — pointing to a systematic rather than random pattern (MNAR).

Median imputation was tried first but immediately rejected. Replacing 276 values with a single median value created a large artificial spike at 0.41 that destroyed the natural shape of the distribution. Random sampling from observed values within each group was used instead, which preserved the original right-skewed shape. A fixed seed (set.seed(123)) was set for reproducibility. Group 0 had 247 NAs filled from 740 observed values, Group 1 had 28 NAs filled from 132 observed values, and the one record with a missing outcome was filled from all observed values.

The before and after histograms confirm the fix worked well. Both show the same right-skewed shape peaking around 0.35 to 0.45, the same range, and no artificial distortion. The mean shifted slightly from 0.4333 to 0.4287 and Q1 and Q3 stayed unchanged at 0.30 and 0.54.

With lipids5 fixed, the full analysis was re-run across all 25 variables producing 50 rows. lipids5 followed lognormal in both groups (Group 0: meanlog = -0.9367, sdlog = 0.4202; Group 1: meanlog = -0.9033, sdlog = 0.4221). Overall results were consistent with the standard task — Group 0 had 19 lognormal and 6 normal, Group 1 had 21 lognormal and 4 normal, and the same 6 variables differed between groups.


