library(tidyverse)

# Question 1: Compare raw and log-transformed protein distributions

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

sample_proteins <- protein_cols[1:4]

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

# Question 2: Identify outliers without trimming

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
  )

outlier_counts %>%
  arrange(desc(outlier_count)) %>%
  print(n = Inf)

group_summary <- outlier_counts %>%
  group_by(group) %>%
  summarise(
    number_of_subjects = n(),
    average_outliers_per_subject = mean(outlier_count),
    median_outliers_per_subject = median(outlier_count),
    subjects_with_outliers = sum(outlier_count > 0),
    .groups = "drop"
  )

print(group_summary)

ggplot(outlier_counts, aes(x = group, y = outlier_count)) +
  geom_boxplot(fill = "gray80", color = "gray30") +
  labs(
    title = "Outlying Protein Measurements per Subject",
    x = "Group",
    y = "Number of Outlying Values per Subject"
  ) +
  theme_minimal()