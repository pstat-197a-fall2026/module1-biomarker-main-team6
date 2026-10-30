library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
load('data/biomarker-clean.RData')

## MULTIPLE TESTING

# Function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- biomarker_clean %>%
  select(-ados) %>%
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  nest(data = c(level, group)) %>% 
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  arrange(p_value) %>%
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# Select significant proteins
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

## RANDOM FOREST

# Store predictors and response separately
predictors <- biomarker_clean %>%
  select(-c(group, ados))

response <- biomarker_clean %>% pull(group) %>% factor()

# Fit RF
set.seed(101422)
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

# Check errors
rf_out$confusion

# Compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

## LOGISTIC REGRESSION

# Select subset of interest
proteins_sstar <- intersect(proteins_s1, proteins_s2)

# Fuzzy intersection:
# using fuzzy matching to grab elements that are similar to each other 
library(stringdist)
proteins_fuzzy <- c()
for (i in 1:length(proteins_s1)) { # compare element-wise similarity of proteins_s1 and proteins_s2
  for (j in 1:length(proteins_s2)) {
    fuzzy_score <- stringdist(proteins_s1[i], proteins_s2[j], method = 'jw')  
    if (fuzzy_score < 0.3) {
      proteins_fuzzy <- c(proteins_fuzzy, proteins_s1[i], proteins_s2[j]) # store both similar elements
    } 
  }
}
proteins_fuzzy <- unique(proteins_fuzzy) # grab the unique proteins
print(proteins_fuzzy)
print(proteins_sstar)

# Discussion: the fuzzy matching algorithm picked up 2 extra proteins, FSTL1 and FSTL3, 
# because they are similar to each other (they are one character apart)

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

biomarker_fuzzy <- biomarker_clean %>%
  select(group, any_of(proteins_fuzzy)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# Partition into training and test set
set.seed(101422)
biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)

biomarker_split_fuzzy <- biomarker_fuzzy %>%
  initial_split(prop = 0.8)

# Fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

fit_fuzzy <- glm(class ~ ., 
                 data = training(biomarker_split_fuzzy), 
                 family = 'binomial')

## RESULTS

# Evaluate errors for both hard intersection + fuzzy intersection on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

# Hard intersection
results <- testing(biomarker_split) %>%
  add_predictions(fit, type = "response") %>%
  mutate(
    actual = factor(class, levels = c(FALSE, TRUE)),
    predicted = factor(pred > 0.5, levels = c(FALSE, TRUE))
  ) 

results %>%
  class_metrics(
    truth = actual,
    estimate = predicted,
    pred,
    event_level = "second"
  )

# Fuzzy intersection
results_fuzzy <- testing(biomarker_split_fuzzy) %>%
  add_predictions(fit_fuzzy, type = "response") %>%
  mutate(
    actual = factor(class, levels = c(FALSE, TRUE)),
    predicted = factor(pred > 0.5, levels = c(FALSE, TRUE))
  ) 

results_fuzzy %>%
  class_metrics(
    truth = actual,
    estimate = predicted,
    pred,
    event_level = "second"
  )

## It seems like the AUC and specificity scores are slightly higher in the fuzzy intersection case 
## compared to the hard intersection case. However, the sensitivity is lower in the fuzzy intersection 
## case than the hard intersection case.

# ROC Curves
library(pROC)
roc_obj <- roc(results$actual, results$pred)

plot(roc_obj, 
     col = "royalblue", 
     main = "ROC Curve (Hard Intersection)", 
     print.auc = TRUE,          
     legacy.axes = TRUE)

roc_obj_fuzzy <- roc(results_fuzzy$actual, results_fuzzy$pred)

plot(roc_obj_fuzzy, 
     col = "royalblue", 
     main = "ROC Curve (Fuzzy Intersection)", 
     print.auc = TRUE,          
     legacy.axes = TRUE)


# Confusion tables 
confusion_table <- table(Predicted = results$predicted, Actual = results$actual)
print(confusion_table)

ggplot(data = as.data.frame(confusion_table), aes(x = Actual, y = Predicted, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  geom_text(aes(label = Freq), size = 6) +
  theme_minimal() +
  labs(title = "Confusion Matrix (Hard Intersection)")

confusion_table_fuzzy <- table(Predicted = results_fuzzy$predicted, Actual = results_fuzzy$actual)
print(confusion_table_fuzzy)

ggplot(data = as.data.frame(confusion_table_fuzzy), aes(x = Actual, y = Predicted, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  geom_text(aes(label = Freq), size = 6) +
  theme_minimal() +
  labs(title = "Confusion Matrix (Fuzzy Intersection)")

## Looking at the confusion matrices, the increase in specificity in the fuzzy intersection case
## is due to an increase in true positives and a decrease in false positives. The decrease in sensitivity,
## on the other hand, is due to a decrease in true positives and an increase in false negatives. 


# Results: While there was an increase in AUC scores when using fuzzy intersection to select the common
# proteins among the RF model and the multiple t-tests, ultimately the magnitude of increase seems very small,
# making it hard to determine if there was any actual statistically significant difference. Additionally, 
# there was a trade-off seen in the logistic regression results with the proteins gathered from fuzzy intersection
# vs. hard intersection: while there was an increase in specificity, there was also a decrease in sensitivity, which 
# makes deciding which of the two sets of proteins to choose more convoluted on top of the small differences. For example, 
# if a higher true positive detection of ASD is the goal, then maybe the proteins chosen by a hard intersection is
# more preferable, etc. So, all in all, it is very hard to tell whether the addition of the 2 proteins made 
# a significant impact to the results before as the differences were small and a clear trade off exists. 









