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
set.seed(1896)
biomarker_split <- biomarker_clean %>%
  initial_split(prop = 0.8, strata = group)

biomarker_train <- training(biomarker_split)
biomarker_test <- testing(biomarker_split)

# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

# whole selection procedure, using the training set only
run_pipeline <- function(n_top){
  
  ## MULTIPLE TESTING
  ####################
  
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
    slice_min(p.adj, n = n_top, with_ties = FALSE) %>%
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
  print(rf_out$confusion)
  
  # compute importance scores
  proteins_s2 <- rf_out$importance %>% 
    as_tibble() %>%
    mutate(protein = rownames(rf_out$importance)) %>%
    slice_max(MeanDecreaseGini, n = n_top, with_ties = FALSE) %>%
    pull(protein)
  
  ## LOGISTIC REGRESSION
  #######################
  
  # select subset of interest
  proteins_sstar <- intersect(proteins_s1, proteins_s2)
  
  if (length(proteins_sstar) == 0) stop('no proteins in common, increase n_top')
  
  prep_fn <- function(.df){
    .df %>%
      select(group, any_of(proteins_sstar)) %>%
      mutate(class = (group == 'ASD')) %>%
      select(-group)
  }
  
  # fit logistic regression model to training set
  fit <- glm(class ~ ., 
             data = prep_fn(biomarker_train), 
             family = 'binomial')
  
  # evaluate errors on test set
  metrics <- prep_fn(biomarker_test) %>%
    add_predictions(fit, type = "response") %>%
    mutate(
      actual = factor(class, levels = c(FALSE, TRUE)),
      predicted = factor(pred > 0.5, levels = c(FALSE, TRUE))
    ) %>%
    class_metrics(
      truth = actual,
      estimate = predicted,
      pred,
      event_level = "second"
    )
  
  list(n_top = n_top, proteins = proteins_sstar, metrics = metrics)
}

res_a <- run_pipeline(n_top = 10)
res_a

# b) choose a larger number (more than ten) of top predictive proteins using each selection method
res_b <- run_pipeline(n_top = 25)
res_b

# compare
map_dfr(list(res_a, res_b), 
        ~ mutate(.x$metrics, n_top = .x$n_top, n_final = length(.x$proteins))) %>%
  select(n_top, n_final, .metric, .estimate) %>%
  pivot_wider(names_from = .metric, values_from = .estimate)