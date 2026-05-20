#--------------------start-------------------------------
# Get current working directory
getwd()
#----------------read dataset--------------------------
install.packages("remotes")
remotes::install_version("wPerm", version = "1.0.1")
library(wPerm)
library(wPerm)
data <- read.csv("C:/Users/DELL/Desktop/Year_1/data analysis/Practical class 4/data_for_analysis.csv")
data$outcome <- as.factor(data$outcome)
summary(data)
# testing for normality of distribution
shapiro.test(data$lipids1)
shapiro.test(data$lipids2)

hist(data$lipids1)  
qqnorm(data$lipids1)

# Spearman's correlation test

spearman_result<-cor.test(data$lipids1, data$lipids2, method="spearman")

print(spearman_result)


# data.frame for result
results <- data.frame(
  variable = character(),
  spearman_corr = numeric(),
  s_p_value = numeric(),
  stringsAsFactors = FALSE
)

# variables for analysis
target_vars <- c("lipids2", "lipids3", "lipids4")

# main 
for (var in target_vars) {
  # spearman
  perm_spearman <- perm.relation(
    x = data$lipids1, 
    y = data[[var]],
    method = "spearman",
    R = 10000
  )
  
  # add result
  results <- rbind(results, data.frame(
    variable = var,
    spearman_corr = perm_spearman$Observed,
    s_p_value =  perm_spearman$p.value
  ))
}


# output result
print(results)

# Task 1 & 2: Spearman correlation with permutation test for hormone variables
hormone_cols <- c("hormone1", "hormone2", "hormone3", "hormone4",
                  "hormone5", "hormone6", "hormone7", "hormone8",
                  "hormone10_generated")

# Spearman correlation matrix
cor_matrix <- cor(data[, hormone_cols], 
                  use = "pairwise.complete.obs", 
                  method = "spearman")

print(cor_matrix)

# Task 2: Permutation test for hormone correlations
hormone_results <- data.frame()

for (i in 1:(length(hormone_cols)-1)) {
  for (j in (i+1):length(hormone_cols)) {
    perm <- perm.relation(
      x = data[[hormone_cols[i]]],
      y = data[[hormone_cols[j]]],
      method = "spearman",
      R = 1000
    )
    hormone_results <- rbind(hormone_results, data.frame(
      Variable1 = hormone_cols[i],
      Variable2 = hormone_cols[j],
      Spearman_rho = round(perm$Observed, 5),
      p_value = perm$p.value
    ))
  }
}


write.csv(hormone_results, "task1_2_hormone_correlation_permutation.csv", row.names = FALSE)


#---Task 3 & 4: Regression analysis between hormone variables---
# Selected pair: hormone3 (predictor) and hormone4 (outcome)
# Reason: hormone3 and hormone4 have the strongest Spearman correlation
# among all hormone pairs (rho = 0.584, p = 0.002), indicating a moderate
# positive relationship. This makes them the most suitable candidates
# for regression analysis among the hormone variables.

df_h <- data
df_h <- df_h[order(df_h$hormone3),]

# linear regression
model_h_linear <- lm(hormone4 ~ hormone3, data = df_h)
summary(model_h_linear)

# second degree polynomial
model_h_2 <- lm(hormone4 ~ poly(hormone3, 2), data = df_h)
summary(model_h_2)

# third degree polynomial
model_h_3 <- lm(hormone4 ~ poly(hormone3, 3), data = df_h)
summary(model_h_3)

# exponential dependence
model_h_exp <- lm(log(hormone4) ~ hormone3, data = df_h)
summary(model_h_exp)

# log dependence
model_h_log <- lm(exp(hormone4) ~ hormone3, data = df_h)
summary(model_h_log)


#---Task 4: BIC model comparison for hormone regression---
rezult_h <- data.frame(
  model = c("model_h_linear", "model_h_2", "model_h_3", 
            "model_h_exp", "model_h_log"),
  BIC_value = c(BIC(model_h_linear), BIC(model_h_2), BIC(model_h_3), 
                BIC(model_h_exp), BIC(model_h_log))
)

rezult_h <- rezult_h[order(rezult_h$BIC_value),]
print(rezult_h)
write.csv(rezult_h, "task4_hormone_BIC_comparison.csv", row.names = FALSE)


#---Task 5: Logistic regression using hormone variables---
# Dependent variable: outcome (binary: 0 = no tumour, 1 = tumour)

# Simple model with one predictor
model_logit_h1 <- glm(outcome ~ hormone1, 
                      data = data, family = binomial)
summary(model_logit_h1)

# Model with two predictors (hormone3 & hormone4 - strongest correlation pair)
model_logit_h2 <- glm(outcome ~ hormone3 + hormone4, 
                      data = data, family = binomial)
summary(model_logit_h2)

# Model with all hormone predictors
model_logit_hall <- glm(outcome ~ hormone1 + hormone2 + hormone3 + 
                          hormone4 + hormone5 + hormone6 + hormone7 + 
                          hormone8 + hormone10_generated, 
                        data = data, family = binomial)
summary(model_logit_hall)

#------visualization of significant results of correlation analysis---------

data<-data[order(data$lipids2),]

plot(data$lipids2, data$lipids1)

lines(data$lipids2, data$lipids1, col = "blue")

abline(lm(data$lipids1 ~ data$lipids2), col="red")




#_____________regression analysis________________ 

df=data
df<-df[order(df$lipids1),]


#linear regression

model_linear <- lm(lipids1 ~ lipids2, data=df)
summary(model_linear)


#second degree polynomal

model_2 <- lm(lipids1 ~ poly(lipids2, 2), data=df)
summary(model_2)

#third degree polynomal

model_3 <- lm(lipids1 ~ poly(lipids2, 3), data=df)

summary(model_3)
#exponential dependence

model_exp <- lm(log(lipids1) ~ lipids2, data=df)
summary(model_exp)
# log dependence

model_log <- lm(exp(lipids1) ~ lipids2, data=df)
summary(model_log)
#comparison of models
#table of result

rezult<-data.frame(model=c("model_linear", "model_2", "model_3", "model_exp", "model_log"), BIC_value=c(BIC(model_linear), BIC(model_2), BIC(model_3), BIC(model_exp), BIC(model_log)))

rezult<-rezult[order(rezult$BIC_value),]

rezult


# __________building graphs______________
#         linear regression graphs

plot(df$lipids2, df$lipids1)
lines(df$lipids2, fitted(model_linear), col="blue")

# Logistic regression
# Dependent variable: outcome 
sum(is.na(data$outcome))  
data <- data[!is.na(data$outcome), ]
# Example: Predicting outcome based on lipids1 and lipids2

# Simple model with one predictor
model_logit_1 <- glm(outcome ~ lipids1, data = data, family = binomial)
summary(model_logit_1)

# Multi-predictor model
model_logit_2 <- glm(outcome ~ lipids1 + lipids2, data = data, family = binomial)
summary(model_logit_2)

# Model with all variables lipids 
model_logit_all <- glm(outcome ~ lipids1 + lipids2 + lipids3 + lipids4, 
                       data = data, family = binomial)
summary(model_logit_all)


# Predicting probabilities for new data (example)
data$pred_prob <- predict(model_logit_2, type = "response")

# Classification by threshold 0.5
data$pred_class <- ifelse(data$pred_prob > 0.5, 1, 0)

# confusion matrix
table(Actual = data$outcome, Predicted = data$pred_class)

# Model Quality Assessment: ROC Curve and AUC 
if (!require(pROC)) install.packages("pROC")
library(pROC)
roc_curve <- roc(data$outcome, data$pred_prob)
plot(roc_curve, main = "ROC-Curve")
auc(roc_curve)

# Stepwise variable selection (AIC)
step_model <- step(model_logit_all, direction = "both")
summary(step_model)
#Coefficient interpretation 
exp(cbind(OR = coef(model_logit_2), confint(model_logit_2)))

#---Task 5 continued: Model evaluation for hormone logistic regression---

# Probability prediction and classification
data$pred_prob_h <- predict(model_logit_hall, type = "response")
data$pred_class_h <- ifelse(data$pred_prob_h > 0.5, 1, 0)

# Confusion matrix
table(Actual = data$outcome, Predicted = data$pred_class_h)

# AIC/BIC comparison of all three hormone models
aic_bic <- data.frame(
  Model = c("model_logit_h1", "model_logit_h2", "model_logit_hall"),
  AIC = c(AIC(model_logit_h1), AIC(model_logit_h2), AIC(model_logit_hall)),
  BIC = c(BIC(model_logit_h1), BIC(model_logit_h2), BIC(model_logit_hall))
)
print(aic_bic)
write.csv(aic_bic, "task5_logistic_AIC_BIC.csv", row.names = FALSE)

# ROC curve and AUC
if (!require(pROC)) install.packages("pROC")
library(pROC)
roc_curve_h <- roc(data$outcome, data$pred_prob_h)
plot(roc_curve_h, main = "ROC Curve - Hormone Logistic Regression")
auc(roc_curve_h)

# Stepwise variable selection
step_model_h <- step(model_logit_hall, direction = "both")
summary(step_model_h)

# Odds ratios and 95% confidence intervals
exp(cbind(OR = coef(model_logit_hall), confint(model_logit_hall)))

# Export confusion matrix
conf_matrix <- as.data.frame(table(Actual = data$outcome, 
                                   Predicted = data$pred_class_h))
write.csv(conf_matrix, "task5_confusion_matrix.csv", row.names = FALSE)

# Export odds ratios
odds_ratios <- as.data.frame(exp(cbind(OR = coef(model_logit_hall), 
                                       confint(model_logit_hall))))
write.csv(odds_ratios, "task5_odds_ratios.csv", row.names = TRUE)

