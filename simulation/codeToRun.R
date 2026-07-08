
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(glmnet)
library(lhs)
library(EmpiricalCalibration)
library(readr)

sampleSize <- c(2, 5) # log
treatmentPrevalence <- c(0.01, 0.99) # lin
numberCovariates <- c(2, 3) # log
unmeasured <- c(-10, -8) # log
outcomePrevalence <- c(-3, -1) # log
treatmentEffect <- c(-1, 1) # lin
numberNCO <- c(0, 3) # log
covariateSd <- c(0, 0.3) # lin
covariatePrev <- c(log(5)/log(10), 3) # log

numberSimulations <- 10
minimumTreatmentPerClassForLasso <- 4

rand <- function(x, range) {
  range[1] + (range[2] - range[1]) * x
}
rand10 <- function(x, range) {
  10 ** (range[1] + (range[2] - range[1]) * x)
}
randPrevalence <- function(n, lambda) {
  u <- runif(n)
  -log(1 - u * (1 - exp(-lambda))) / lambda
}

parameters <- randomLHS(numberSimulations, 9)

parameters <- tibble(
  sample_size = round(rand10(parameters[,1], sampleSize)),
  treatment_prevalence = rand(parameters[,2], treatmentPrevalence),
  number_covariates = round(rand10(parameters[,3], numberCovariates)),
  unmeasured = round(rand10(parameters[,4], unmeasured) * number_covariates),
  outcome_prevalence = rand10(parameters[,5], outcomePrevalence),
  treatment_effect = rand(parameters[,6], treatmentEffect),
  nco = round(rand10(parameters[,7], numberNCO)),
  covariate_sd = rand(parameters[,8], covariateSd),
  covariate_prev = rand10(parameters[,9], covariatePrev)
) |>
  filter(.data$sample_size > 2 * .data$number_covariates)

fmt <- paste0("cov_%0", ceiling(numberCovariates[2]),"i")
fmtnco <- paste0("nco_%0", ceiling(numberCovariates[2]),"i")

result <- parameters |>
  pmap(\(sample_size, treatment_prevalence, number_covariates, unmeasured, outcome_prevalence, treatment_effect, nco, covariate_sd, covariate_prev) {
    ts <- Sys.time()
    n <- sample_size

    # result
    result <- tibble(
      sample_size = n,
      treatment_prevalence = treatment_prevalence,
      number_covariates = number_covariates,
      unmeasured = unmeasured,
      outcome_prevalence = outcome_prevalence,
      treatment_effect = treatment_effect,
      nco = nco,
      covariate_sd = covariate_sd,
      covariate_prev = 1/covariate_prev
    )

    # covariates
    covariates <- randPrevalence(n = number_covariates, lambda = covariate_prev) |>
      imap(\(p, i) {
        tibble(!!sprintf(fmt, i) := if_else(runif(n = n) < p, 1, 0))
      }) |>
      bind_cols()

    # covariates -> treatment
    cov <- rnorm(n = number_covariates, sd = covariate_sd)

    # treatment
    p_treatment <- as.vector(as.matrix(covariates) %*% cov)
    fn <- \(x) mean(1 / (1 + exp(-(x + p_treatment)))) - treatment_prevalence
    delta <- uniroot(fn, interval = c(-100, 100))$root
    p_treatment <- 1 / (1 + exp(-(delta + p_treatment)))
    treatment <- if_else(runif(n = n) < p_treatment, 1, 0)

    treatment_class_count <- table(treatment)
    if (length(treatment_class_count) != 2) {
      return(NULL)
    }

    # covariates -> outcome
    cov <- rnorm(n = number_covariates, sd = covariate_sd)

    # outcome
    linear_outcome <- as.vector(as.matrix(covariates) %*% cov)
    linear_outcome_observed <- linear_outcome + treatment_effect * treatment
    fn <- \(x) mean(1 / (1 + exp(-(x + linear_outcome_observed)))) - outcome_prevalence
    delta <- uniroot(fn, interval = c(-100, 100))$root
    p_outcome_0 <- 1 / (1 + exp(-(delta + linear_outcome)))
    p_outcome_1 <- 1 / (1 + exp(-(delta + linear_outcome + treatment_effect)))
    p_outcome <- if_else(treatment == 1, p_outcome_1, p_outcome_0)
    outcome <- if_else(runif(n = n) < p_outcome, 1, 0)

    marginal_treatment_effect <- qlogis(mean(p_outcome_1)) - qlogis(mean(p_outcome_0))
    overlap_target_weight <- p_treatment * (1 - p_treatment)
    overlap_treatment_effect <- qlogis(weighted.mean(p_outcome_1, overlap_target_weight)) -
      qlogis(weighted.mean(p_outcome_0, overlap_target_weight))
    result$marginal_treatment_effect <- marginal_treatment_effect
    result$overlap_treatment_effect <- overlap_treatment_effect

    # negative control outcomes
    p_nco <- as.vector(as.matrix(covariates) %*% cov)
    nco <- 10 ** runif(n = nco, min = -3, max = -0.5) |>
      imap(\(p, i) {
        fn <- \(x) mean(1 / (1 + exp(-(x + p_nco)))) - p
        delta <- uniroot(fn, interval = c(-100, 100))$root
        p_nco <- 1 / (1 + exp(-(delta + p_nco)))
        tibble(!!sprintf(fmtnco, i) := if_else(runif(n = n) < p_nco, 1, 0))
      }) |>
      bind_cols()

    # remove unmeasured confounding
    covariates <- covariates |>
      select(!seq_len(unmeasured))

    # crude coefficient
    x <- tibble(outcome = outcome, treatment = treatment)
    crude <- glm(outcome ~ treatment, data = x, family = binomial(link = "logit")) |>
      tidy() |>
      filter(term == "treatment")
    result$crude_error <- crude$estimate - marginal_treatment_effect
    result$crude_sd <- crude$std.error

    # lasso regression
    if (ncol(covariates) > 0 && min(treatment_class_count) >= minimumTreatmentPerClassForLasso) {
      cv_fit <- cv.glmnet(
        x = as.matrix(covariates),
        y = treatment,
        family = "binomial",
        alpha = 1,
        nfolds = min(10, min(treatment_class_count))
      )
      lasso_coef <- coef(cv_fit, s = "lambda.min")
      covs <- rownames(lasso_coef)[as.vector(lasso_coef != 0)]
      covs <- covs[startsWith(covs, "cov_")]
    } else {
      covs <- character()
    }

    # weights
    if (length(covs) > 0) {
      formula <- paste0("treatment ~ ", paste0(covs, collapse = " + "))
      fit <- glm(formula = formula, family = binomial(link = "logit"), data = bind_cols(x, covariates))
      prob <- predict(fit, type = "response")
    } else {
      prob <- mean(treatment)
    }
    
    x <- x |>
      mutate(
        ps = prob,
        weight = if_else(treatment == 1, 1 - ps, ps)
      )

    # weighted coefficient
    weighted <- glm(outcome ~ treatment, data = x, weights = weight, family = binomial(link = "logit")) |>
      tidy() |>
      filter(term == "treatment")
    result$weighted_error <- weighted$estimate - overlap_treatment_effect
    result$weighted_sd <- weighted$std.error

    # negative control outcomes
    ncos <- nco |>
      map(\(out) {
        x$outcome <- out
        glm(outcome ~ treatment, data = x, weights = weight, family = binomial(link = "logit")) |>
          tidy() |>
          filter(term == "treatment") |>
          select("estimate", "std" = "std.error")
      }) |>
      bind_rows()
    model <- fitSystematicErrorModel(ncos$estimate, ncos$std, rep(0, nrow(ncos)))
    result$nco_model_mean_intercept <- model[1]
    result$nco_model_mean_slope <- model[2]
    result$nco_model_sd_intercept <- model[3]
    result$nco_model_sd_slope <- model[4]
    calibrated <- calibrateConfidenceInterval(weighted$estimate, weighted$std.error, model)
    result$calibrated_error <- calibrated$logRr - overlap_treatment_effect
    result$calibrated_sd <- calibrated$seLogRr

    # adjusted model
    x <- x |>
      bind_cols(nco)
    formula <- paste0("outcome ~ treatment + ", paste0(colnames(nco), collapse = " + "))
    crude_adjusted <- glm(formula = formula, data = x, family = binomial(link = "logit")) |>
      tidy() |>
      filter(term == "treatment")
    result$crude_adjusted_error <- crude_adjusted$estimate - marginal_treatment_effect
    result$crude_adjusted_sd <- crude_adjusted$std.error
    weighted_adjusted <- glm(formula = formula, data = x, weights = weight, family = binomial(link = "logit")) |>
      tidy() |>
      filter(term == "treatment")
    result$weighted_adjusted_error <- weighted_adjusted$estimate - overlap_treatment_effect
    result$weighted_adjusted_sd <- weighted_adjusted$std.error

    # measured
    result$measured_treatment_prevalence <- mean(x$treatment)
    result$measured_outcome_prevalence <- mean(x$outcome)

    # time
    result$time <- as.numeric(difftime(time1 = Sys.time(), time2 = ts))

    return(result)
  }) |>
  bind_rows()

write_csv(x = result, file = "simulations.csv")
