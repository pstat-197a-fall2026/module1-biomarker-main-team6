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
