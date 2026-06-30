
library(cli)
library(glmnet)
library(dplyr)
library(purrr)
library(rlang)
library(readr)
library(tidyr)
library(here)
library(EmpiricalCalibration)

# load data
data <- read_csv(file = here("clean_data.csv"), show_col_types = FALSE) |>
  mutate(
    hydroxychloroquine = if_else(exposure == "hydroxychloroquine", 1, 0),
    sulfasalazine = if_else(exposure == "sulfasalazine", 1, 0),
    leflunomide = if_else(exposure == "leflunomide", 1, 0)
  )

# source functions
source("functions.R")

# Define parameters
number_simulations <- 50

number_confounders <- c(50)
unmeasured_fraction <- c(0, 0.25, 0.5, 0.75, 1)
sample_size <- c(5000, 10000, 20000)
outcome_prev <- c(0.2)
negative_control_outcomes <- 50
treatment_effect_hydroxychloroquine <- 0.3
treatment_effect_sulfasalazine <- 0.2
treatment_effect_leflunomide <- 0.1

# expand parameters
params <- expand_grid(
  number_confounders = number_confounders,
  unmeasured_fraction = unmeasured_fraction,
  sample_size = sample_size,
  outcome_prev = outcome_prev,
  number_simulations = number_simulations,
  negative_control_outcomes = negative_control_outcomes,
  treatment_effect_hydroxychloroquine = treatment_effect_hydroxychloroquine,
  treatment_effect_sulfasalazine = treatment_effect_sulfasalazine,
  treatment_effect_leflunomide = treatment_effect_leflunomide
)

# prepare outcome models
outcome_models <- list()
for (nc in number_confounders) {
  outcome_models[[sprintf("%i", nc)]] <- fit_outcome_model(data = data, number_confounders = nc)
}

# lasso for selected variables
variables_lasso <- list()
for (nm in c("hydroxychloroquine", "sulfasalazine", "leflunomide")) {
  data_lasso <- data |>
    filter(exposure %in% c("methotrexate", nm)) |>
    mutate(exposure = as.numeric(exposure == nm)) |>
    select("exposure", starts_with("var_"))
  lasso <- cv.glmnet(
    x = as.matrix(select(data_lasso, starts_with("var_"))),
    y = data_lasso$exposure,
    alpha = 1,
  )
  best_lambda <- lasso$lambda.min
  coefs <- coef(lasso, s = "lambda.min")
  variables_lasso[[nm]] <- as.numeric(coefs) |>
    set_names(rownames(coefs)) |>
    keep(\(x) x != 0) |>
    names() |>
    keep(\(x) startsWith(x = x, prefix = "var_"))
}

# subset data
data <- data |>
  select("id", "exposure", "age", "sex", any_of(unique(unlist(variables_lasso))))

result <- list()
for (k in seq_len(nrow(params))) {
  # get outcome models and confounders_candidates
  nc <- sprintf("%i", params$number_confounders[k])
  confounders_candidates <- outcome_models[[nc]]$confounders_candidates
  outcome_model <- outcome_models[[nc]]$outcome_model

  # treatment effect
  treatment_effect <- list(
    hydroxychloroquine = params$treatment_effect_hydroxychloroquine[k],
    sulfasalazine = params$treatment_effect_sulfasalazine[k],
    leflunomide = params$treatment_effect_leflunomide[k]
  )

  cli_inform("Parameter combination: {k} out of {nrow(params)}")
  result[[k]] <- compute_plasmode_data_real_exposure(
    data = data,
    treatment_effect = treatment_effect,
    outcome_model = outcome_model,
    confounders_candidates = confounders_candidates,
    variables_lasso = variables_lasso,
    sample_size = params$sample_size[k],
    outcome_prev = params$outcome_prev[k],
    number_simulations = params$number_simulations[k],
    negative_control_outcomes = params$negative_control_outcomes[k],
    unmeasured_fraction = params$unmeasured_fraction[k]
  )
}
result <- bind_rows(result)

x <- result |>
  mutate(error = abs(coef_real - coef)) |>
  group_by(type, comparator, unmeasured) |>
  summarise(
    y = median(error),
    y_max = quantile(error, 0.75),
    y_min = quantile(error, 0.25),
    .groups = "drop"
  )

library(ggplot2)

ggplot(data = x, mapping = aes(x = unmeasured, y = y, ymax = y_max, ymin = y_min, colour = type)) +
  geom_point() +
  geom_errorbar() +
  facet_grid(comparator ~ .)
