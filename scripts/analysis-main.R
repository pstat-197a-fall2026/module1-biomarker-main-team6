library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
library(stringdist)
load('data/biomarker-clean.RData')
# Question 1

var_names <- read_csv(
  "data/biomarker-raw.csv",
  col_names = FALSE,
  n_max = 2,
  col_select = -(1:2),
  show_col_types = FALSE
) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, abbreviation = V2) %>%
  na.omit()

raw <- read_csv(
  "data/biomarker-raw.csv",
  skip = 2,
  col_select = -2L,
  col_names = c(
    "group",
    "empty",
    pull(var_names, abbreviation),
    "ados"
  ),
  na = c("-", ""),
  show_col_types = FALSE
) %>%
  filter(!is.na(group))

protein_cols <- setdiff(names(raw), c("group", "empty", "ados"))

set.seed(1)
sample_proteins <- sample(protein_cols, 4)

for (protein in sample_proteins) {
  values <- raw[[protein]]
  values <- values[is.finite(values) & values > 0]
  
  hist(
    values,
    main = paste("Raw Protein Levels:", protein),
    xlab = "Protein Level",
    ylab = "Number of Subjects",
    col = "gray80",
    border = "white",
    breaks = 20
  )
  
  hist(
    log10(values),
    main = paste("Log-Transformed Protein Levels:", protein),
    xlab = "Log10(Protein Level)",
    ylab = "Number of Subjects",
    col = "gray80",
    border = "white",
    breaks = 20
  )
}

# Question 2

biomarker_untrimmed <- raw %>%
  mutate(
    subject = row_number(),
    across(
      all_of(protein_cols),
      ~ ifelse(.x > 0, log10(.x), NA_real_)
    )
  ) %>%
  mutate(
    across(
      all_of(protein_cols),
      ~ as.numeric(scale(.x))
    )
  ) %>%
  select(subject, group, all_of(protein_cols))

outlier_counts <- biomarker_untrimmed %>%
  pivot_longer(
    cols = all_of(protein_cols),
    names_to = "protein",
    values_to = "value"
  ) %>%
  group_by(subject, group) %>%
  summarise(
    total_measurements = sum(!is.na(value)),
    outlier_count = sum(abs(value) > 3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    outlier_rate = outlier_count / total_measurements
  )

outlier_counts %>%
  arrange(desc(outlier_count)) %>%
  head(15) %>%
  print()

group_summary <- outlier_counts %>%
  group_by(group) %>%
  summarise(
    number_of_subjects = n(),
    average_outliers_per_subject = mean(outlier_count),
    median_outliers_per_subject = median(outlier_count),
    subjects_with_outliers = sum(outlier_count > 0),
    average_outlier_rate = mean(outlier_rate),
    .groups = "drop"
  )

print(group_summary)

wilcox.test(outlier_count ~ group, data = outlier_counts)

ggplot(outlier_counts, aes(x = group, y = outlier_count)) +
  geom_boxplot(fill = "gray80", color = "gray30", outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
  labs(
    title = "Outlying Protein Measurements per Subject",
    x = "Group",
    y = "Number of Outlying Values per Subject"
  ) +
  theme_minimal()

# Question 3
library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
load('data/biomarker-clean.RData')

library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
load('data/biomarker-clean.RData')

# a) repeat the analysis but carry out the entire selection procedure on a training partition 
# in other words, set aside some testing data at the very beginning and don't use 
# it until you are evaluating accuracy at the very end

# I will choose to use a 80/20 split as that is the most common split

# partition into training and test set
set.seed(1896)
biomarker_split <- biomarker_clean %>%
  initial_split(prop = 0.8)

biomarker_train <- training(biomarker_split)
biomarker_test <- testing(biomarker_split)

## MULTIPLE TESTING
####################

# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- biomarker_train %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# select significant proteins
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

## RANDOM FOREST
##################

# store predictors and response separately
predictors <- biomarker_train %>%
  select(-c(group, ados))

response <- biomarker_train %>% pull(group) %>% factor()

# fit RF
set.seed(101422)
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

# check errors
rf_out$confusion

# compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

## LOGISTIC REGRESSION
#######################

# select subset of interest
proteins_sstar <- intersect(proteins_s1, proteins_s2)

biomarker_sstar_train <- biomarker_train %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

biomarker_sstar_test <- biomarker_test %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = biomarker_sstar_train, 
           family = 'binomial')

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

preds_a <- biomarker_sstar_test %>%
  add_predictions(fit, type = "response") %>%
  mutate(
    actual = factor(class, levels = c(FALSE, TRUE)),
    predicted = factor(pred > 0.5, levels = c(FALSE, TRUE))
  )

results_a <- preds_a %>%
  class_metrics(
    truth = actual,
    estimate = predicted,
    pred,
    event_level = "second"
  )

proteins_sstar_a <- proteins_sstar
proteins_sstar_a
results_a


# b) choose a larger number (more than ten) of top predictive proteins using each selection method
# same as above, only n = 10 changed to n = 25 
# (t tests and RF are already fit on the training set, so only the selection changes)

## MULTIPLE TESTING
####################

# select significant proteins
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 25) %>%
  pull(protein)

## RANDOM FOREST
##################

# compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 25) %>%
  pull(protein)

## LOGISTIC REGRESSION
#######################

# select subset of interest
proteins_sstar <- intersect(proteins_s1, proteins_s2)

biomarker_sstar_train <- biomarker_train %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

biomarker_sstar_test <- biomarker_test %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = biomarker_sstar_train, 
           family = 'binomial')

# evaluate errors on test set
preds_b <- biomarker_sstar_test %>%
  add_predictions(fit, type = "response") %>%
  mutate(
    actual = factor(class, levels = c(FALSE, TRUE)),
    predicted = factor(pred > 0.5, levels = c(FALSE, TRUE))
  )

results_b <- preds_b %>%
  class_metrics(
    truth = actual,
    estimate = predicted,
    pred,
    event_level = "second"
  )

proteins_sstar_b <- proteins_sstar
proteins_sstar_b
results_b


## COMPARE RESULTS
###################

# label each run with the number of proteins in the final model
run_labels <- c(
  paste0('Top 10 per method (', length(proteins_sstar_a), ' in final model)'),
  paste0('Top 25 per method (', length(proteins_sstar_b), ' in final model)')
)

preds_all <- bind_rows(preds_a, preds_b, .id = 'run') %>%
  mutate(run = factor(run, levels = c('1', '2'), labels = run_labels))

# table of test set metrics
preds_all %>%
  group_by(run) %>%
  class_metrics(
    truth = actual,
    estimate = predicted,
    pred,
    event_level = "second"
  ) %>%
  select(run, .metric, .estimate) %>%
  mutate(.estimate = round(.estimate, 3)) %>%
  pivot_wider(names_from = .metric, values_from = .estimate)

# ROC curves
auc_labels <- preds_all %>%
  group_by(run) %>%
  roc_auc(truth = actual, pred, event_level = 'second') %>%
  mutate(label = paste0('AUC: ', round(.estimate, 3)))

preds_all %>%
  group_by(run) %>%
  roc_curve(truth = actual, pred, event_level = 'second') %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_abline(color = 'gray') +
  geom_path(color = 'royalblue', linewidth = 0.8) +
  geom_text(data = auc_labels, 
            aes(x = 0.65, y = 0.35, label = label), 
            color = 'royalblue', 
            inherit.aes = FALSE) +
  facet_wrap(~ run) +
  coord_equal() +
  labs(title = 'ROC Curve', 
       x = '1 - Specificity', 
       y = 'Sensitivity') +
  theme_bw()

# c)
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

# Question 4
# Task 4: smaller panels versus the in-class model.
# This comparison inherits the in-class full-data preprocessing and screening.
source("scripts/inclass-analysis.R")

train_data <- training(biomarker_split)
test_data <- testing(biomarker_split)

# Fit logistic regression and return predicted ASD probabilities.
predict_panel <- function(train, test, proteins) {
  model <- glm(class ~ ., family = binomial(),
               data = train %>% select(class, all_of(proteins)))
  tibble(
    actual = factor(test$class, levels = c(FALSE, TRUE)),
    pred = as.numeric(predict(model, newdata = test, type = "response"))
  )
}

# Try all 31 nonempty subsets of the five in-class proteins.
panels <- unlist(lapply(seq_along(proteins_sstar), function(size) {
  combn(proteins_sstar, size, simplify = FALSE)
}), recursive = FALSE)

# Use the same stratified five-fold splits for every panel, repeated 10 times.
set.seed(20261029)
folds <- vfold_cv(
  train_data %>% mutate(class = factor(class, levels = c(FALSE, TRUE))),
  v = 5, repeats = 10, strata = class
)

cv_results <- map_dfr(panels, function(protein_set) {
  predictions <- map2_dfr(folds$splits, folds$id, function(split, repeat_id) {
    predict_panel(analysis(split), assessment(split), protein_set) %>%
      mutate(repeat_id = repeat_id)
  })
  
  predictions %>%
    group_by(repeat_id) %>%
    summarise(
      accuracy = mean((pred > 0.5) == (actual == "TRUE")),
      auc = roc_auc_vec(actual, pred, event_level = "second"),
      .groups = "drop"
    ) %>%
    summarise(cv_accuracy = mean(accuracy), cv_auc = mean(auc)) %>%
    mutate(proteins = list(protein_set),
           n_proteins = length(protein_set),
           panel = paste(protein_set, collapse = ", "))
})

# Select candidates using training CV, before evaluating their test results.
baseline <- cv_results %>% filter(n_proteins == length(proteins_sstar))

# "Comparable" means at most 0.02 lower CV accuracy and AUC than baseline.
smallest <- cv_results %>%
  filter(cv_accuracy >= baseline$cv_accuracy - 0.02,
         cv_auc >= baseline$cv_auc - 0.02) %>%
  arrange(n_proteins, desc(cv_accuracy), desc(cv_auc)) %>%
  slice_head(n = 1)

best <- cv_results %>%
  arrange(desc(cv_accuracy), desc(cv_auc), n_proteins) %>%
  slice_head(n = 1)

selected <- bind_rows(
  baseline %>% mutate(model = "In-class baseline"),
  smallest %>% mutate(model = "Smallest CV-comparable panel"),
  best %>% mutate(model = "Best training-CV panel")
)

# Compare the chosen panels on the same test subjects at a 0.5 threshold.
comparison <- selected %>%
  mutate(metrics = map(proteins, function(protein_set) {
    predict_panel(train_data, test_data, protein_set) %>%
      mutate(predicted = factor(pred > 0.5, levels = c(FALSE, TRUE))) %>%
      class_metrics(truth = actual, estimate = predicted, pred,
                    event_level = "second") %>%
      select(.metric, .estimate) %>%
      pivot_wider(names_from = .metric, values_from = .estimate)
  })) %>%
  unnest(metrics) %>%
  select(model, panel, n_proteins, cv_accuracy, cv_auc,
         accuracy, sensitivity, specificity, roc_auc)

# Report tables: all performance measures are proportions on a 0-1 scale.
cv_results %>%
  group_by(n_proteins) %>%
  arrange(desc(cv_accuracy), desc(cv_auc), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(Proteins = n_proteins, Panel = panel,
         `CV accuracy` = cv_accuracy, `CV AUC` = cv_auc) %>%
  knitr::kable(format = "pipe", digits = 3,
               caption = "Best training cross-validation accuracy at each panel size.") %>%
  print()

comparison %>%
  select(Panel = panel, Proteins = n_proteins, Accuracy = accuracy,
         Sensitivity = sensitivity, Specificity = specificity, AUC = roc_auc) %>%
  knitr::kable(format = "pipe", digits = 3,
               caption = "Test-set performance on the same 31 subjects (threshold = 0.5).") %>%
  print()
