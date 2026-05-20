#--------------------start-------------------------------
# Get current working directory
getwd()
#----------------read dataset--------------------------
data_for_analysis<-read.csv("C:/Users/DELL/Desktop/Year_1/data analysis/Practical class 3/data_for_analysis.csv")

#-----------descriptive statistics------------------
summary(data_for_analysis)

#-----------for publication tables-----------------
#---------------Creating a custom table--------------
# Homework: Creating a custom table with descriptive statistics results

install.packages("gtsummary")
install.packages(c("cardx", "cards"))
library(cardx)
library(gtsummary)

tbl_summary(data_for_analysis)  # Automatic table
tbl_summary(data_for_analysis, by = outcome)  # By groups

#--------------Statistical Tests---------------------
value_outcome1<-data_for_analysis[data_for_analysis$outcome=="1",]$lipids1
hist(value_outcome1, col = "lightblue")

qqnorm(value_outcome1, main = "Q-Q Plot")
qqline(value_outcome1, col = "red", lwd = 2)

# Shapiro-Wilk test (for n < 5000)
shapiro.test(value_outcome1)


value_outcome0<-data_for_analysis[data_for_analysis$outcome=="0",]$lipids1
hist(value_outcome0, col = "lightgreen")

qqnorm(value_outcome0, main = "Q-Q Plot")
qqline(value_outcome0, col = "red", lwd = 2)

# Shapiro-Wilk test (for n < 5000)
shapiro.test(value_outcome0)

#-------Levene's Test for Homogeneity of Variance--------------
install.packages("car")
library(car)
str(data_for_analysis)
data_for_analysis$outcome<- as.factor(data_for_analysis$outcome)
car::leveneTest(lipids1 ~ outcome, data = data_for_analysis)
#---------------Application of the Brunner-Munzel test----------
install.packages("lawstat")
library(lawstat)
group1 <- data_for_analysis$lipids1[data_for_analysis$outcome == "0"]
group2 <- data_for_analysis$lipids1[data_for_analysis$outcome == "1"]

brunner.munzel.test(group1, group2)
#-------------comparison of results with other tests--------------
t.test(group1, group2)
wilcox.test(group1, group2)

#-------------comparison of results with other tests--------------
t.test(group1, group2)
wilcox.test(group1, group2)

#-------------------define hormone columns------------------------
hormone_cols <- c("hormone1", "hormone2", "hormone3", "hormone4",
                  "hormone5", "hormone6", "hormone7", "hormone8",
                  "hormone10_generated")


#-------------------Task 1: Hormone descriptive table with p-values----
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

#-------------------Task 2: Shapiro-Wilk + Levene's Test----------
sw_results <- data.frame()

for (h in hormone_cols) {
  for (grp in c("0", "1")) {
    vals <- data_for_analysis[[h]][data_for_analysis$outcome == grp]
    vals <- vals[!is.na(vals)]
    sw <- shapiro.test(vals)
    sw_results <- rbind(sw_results, data.frame(
      Hormone = h,
      Group   = grp,
      W       = round(sw$statistic, 5),
      p_value = round(sw$p.value, 5),
      Normal  = ifelse(sw$p.value >= 0.05, "Yes", "No"),
      stringsAsFactors = FALSE
    ))
  }
}
write.csv(sw_results, "shapiro_wilk_test_results.csv", row.names = FALSE)


#-------------------Levene's Test for all hormones----------------
lev_results <- data.frame()

for (h in hormone_cols) {
  lev <- car::leveneTest(data_for_analysis[[h]] ~ data_for_analysis$outcome,
                         center = median)
  lev_results <- rbind(lev_results, data.frame(
    Hormone   = h,
    F_value   = round(lev$`F value`[1], 4),
    p_value   = round(lev$`Pr(>F)`[1], 5),
    Equal_Var = ifelse(lev$`Pr(>F)`[1] >= 0.05, "Yes", "No"),
    stringsAsFactors = FALSE
  ))
}

write.csv(lev_results, "levene_test_results.csv", row.names = FALSE)


#-------------------Task 3: Histograms and Q-Q Plots for all hormones-----------
#-------------------Histograms for all hormones-----------
for (h in hormone_cols) {
  par(mfrow = c(1, 2))
  
  g0 <- data_for_analysis[[h]][data_for_analysis$outcome == "0"]
  g1 <- data_for_analysis[[h]][data_for_analysis$outcome == "1"]
  
  hist(g0, col = "lightblue", main = paste(h, "| Outcome = 0"),
       xlab = h, breaks = 30)
  hist(g1, col = "salmon", main = paste(h, "| Outcome = 1"),
       xlab = h, breaks = 30)
}

par(mfrow = c(1, 1))


#-------------------Q-Q Plots for all hormones-----------
for (h in hormone_cols) {
  par(mfrow = c(1, 2))
  
  g0 <- data_for_analysis[[h]][data_for_analysis$outcome == "0"]
  g1 <- data_for_analysis[[h]][data_for_analysis$outcome == "1"]
  
  qqnorm(g0, main = paste("Q-Q:", h, "| Outcome = 0"), col = "steelblue")
  qqline(g0, col = "red", lwd = 2)
  
  qqnorm(g1, main = paste("Q-Q:", h, "| Outcome = 1"), col = "darkorange")
  qqline(g1, col = "red", lwd = 2)
}

par(mfrow = c(1, 1))



#---Task 4: Brunner-Munzel, t.test, Wilcox for 2 independent groups (all hormones)----

test_results <- data.frame()

for (h in hormone_cols) {
  g0 <- data_for_analysis[[h]][data_for_analysis$outcome == "0"]
  g1 <- data_for_analysis[[h]][data_for_analysis$outcome == "1"]
  
  g0 <- g0[!is.na(g0)]
  g1 <- g1[!is.na(g1)]
  
  bm <- brunner.munzel.test(g0, g1)
  tt <- t.test(g0, g1)
  wt <- wilcox.test(g0, g1)
  
  test_results <- rbind(test_results, data.frame(
    Hormone          = h,
    Median_IQR_G0    = paste0(round(median(g0), 3), " (", round(quantile(g0, 0.25), 3), ", ", round(quantile(g0, 0.75), 3), ")"),
    Median_IQR_G1    = paste0(round(median(g1), 3), " (", round(quantile(g1, 0.25), 3), ", ", round(quantile(g1, 0.75), 3), ")"),
    p_BrunnerMunzel  = round(bm$p.value, 5),
    p_ttest          = round(tt$p.value, 5),
    p_Wilcoxon       = round(wt$p.value, 5),
    stringsAsFactors = FALSE
  ))
}

print(test_results)
write.csv(test_results, "task4_descriptive_with_pvalues.csv", row.names = FALSE)


#Conclusion: which test is applicable?

#The Wilcoxon rank-sum test is applicable.
#Because the data is non-normal (Shapiro-Wilk, p < 0.05) and variances are equal between groups (Levene's test, p > 0.05).

#---Task 5: Correlation heatmap for all hormones by group (Spearman)----
install.packages(c("ggplot2", "reshape2"))
library(ggplot2)
library(reshape2)

# Group 0 - no tumour
df_g0 <- data_for_analysis[data_for_analysis$outcome == "0", hormone_cols]
cor_g0 <- cor(df_g0, use = "pairwise.complete.obs", method = "spearman")

# Group 1 - tumour
df_g1 <- data_for_analysis[data_for_analysis$outcome == "1", hormone_cols]
cor_g1 <- cor(df_g1, use = "pairwise.complete.obs", method = "spearman")

# Plot Group 0
melted_g0 <- melt(cor_g0)
p_g0 <- ggplot(melted_g0, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 2)), size = 3) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman r") +
  labs(title = "Hormone Correlation Heatmap - Outcome 0 (No Tumour)",
       x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot Group 1
melted_g1 <- melt(cor_g1)
p_g1 <- ggplot(melted_g1, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 2)), size = 3) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman r") +
  labs(title = "Hormone Correlation Heatmap - Outcome 1 (Tumour)",
       x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print plots
# Save plots as PNG
ggsave("task5_correlation_heatmap_outcome0.png", plot = p_g0, 
       width = 8, height = 6, dpi = 300)
ggsave("task5_correlation_heatmap_outcome1.png", plot = p_g1, 
       width = 8, height = 6, dpi = 300)

# Export correlation matrices as CSV
write.csv(cor_g0, "task5_correlation_matrix_outcome0.csv", row.names = TRUE)
write.csv(cor_g1, "task5_correlation_matrix_outcome1.csv", row.names = TRUE)


#----------------------------EDA----------------------------------
install.packages("DataExplorer")
library(DataExplorer)
create_report(data_for_analysis)  # Generates HTML report with graphs and statistics
create_report(
  data = data_for_analysis,
  output_file = "EDA_Report.html",  
  output_dir = getwd(),                
  report_title = "EDA Report"          
)

names(data_for_analysis)

