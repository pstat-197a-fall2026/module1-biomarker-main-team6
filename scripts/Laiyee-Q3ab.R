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
