
compute_plasmode_data_real_exposure <- function(data,
                                                treatment_effect,
                                                outcome_model,
                                                confounders_candidates,
                                                variables_lasso,
                                                sample_size,
                                                outcome_prev,
                                                number_simulations,
                                                negative_control_outcomes,
                                                unmeasured_fraction) {
  cli_inform(c(i = "Prepare parameters."))
  time <- list(lasso = 0, weights = 0, crude = 0, calibrated = 0, adjusted = 0)
  ti <- as.numeric(Sys.time())

  # adjust for desire treatment effect
  x <- model.matrix(outcome_model)
  coeff <- coef(outcome_model)
  for (nm in names(treatment_effect)) {
    coeff[[nm]] <- treatment_effect[[nm]]
  }
  x_outcome <- as.vector(x %*% coeff)
  fn <- \(x) mean(1 / (1 + exp(-(x + x_outcome)))) - outcome_prev
  delta <- uniroot(fn, interval = c(-20, 20))$root
  outcome_prob <- 1 / (1 + exp(-(delta + x_outcome)))

  # adjust the negative control treatment effect
  coeff <- coef(outcome_model)
  coeff[names(coeff) %in% names(treatment_effect)] <- 0
  x_nco <- as.vector(x %*% coeff)

  # result
  result <- list()

  # simulate data
  for (i in seq_len(number_simulations)) {
    cli_inform(c("simulation: {i} out of {number_simulations}"))

    # sample data
    id <- sample(x = seq_len(nrow(data)), size = sample_size, replace = FALSE)
    sub_data <- data[id,]

    # outcome
    out <- list()
    out[[1]] <- rbinom(n = sample_size, size = 1, prob = outcome_prob[id])

    # negative control outcomes
    for (j in seq_len(negative_control_outcomes)) {
      nco_prev <- rbeta(n = 1, shape1 = 3, shape2 = 10)
      fn <- \(x) mean(1 / (1 + exp(-(x + x_nco)))) - nco_prev
      delta <- uniroot(fn, interval = c(-20, 20))$root
      out[[j + 1]] <- rbinom(n = sample_size, size = 1, prob = 1 / (1 + exp(-(delta + x_nco[id]))))
    }
    names(out) <- c("outcome", sprintf("nco_%03i", seq_len(negative_control_outcomes)))

    # dataset
    sub_data <- sub_data |>
      mutate(!!!out)

    # unmeasured confounding
    unmeasured <- confounders_candidates[runif(n = length(confounders_candidates)) <= unmeasured_fraction]
    sub_data <- sub_data |>
      select(!any_of(unmeasured))

    # loop for exposures
    result_i <- list()
    for (nm in c("hydroxychloroquine", "sulfasalazine", "leflunomide")) {
      result_nm <- list()

      # filter to exposures of interest
      x <- sub_data |>
        filter(exposure %in% c("methotrexate", nm)) |>
        mutate(exposure = as.numeric(exposure == nm)) |>
        select("exposure", "outcome", starts_with("nco_"), any_of(c("age", "sex", variables_lasso[[nm]])))

      # get variables
      t0 <- as.numeric(Sys.time())
      lasso <- cv.glmnet(
        x = as.matrix(select(x, starts_with("var_"))),
        y = x$exposure,
        alpha = 1,
      )
      best_lambda <- lasso$lambda.min
      coefs <- coef(lasso, s = "lambda.min")
      variables <- as.numeric(coefs) |>
        set_names(rownames(coefs)) |>
        keep(\(x) x != 0) |>
        names() |>
        keep(\(x) startsWith(x = x, prefix = "var_"))
      time$lasso <- time$lasso + as.numeric(Sys.time()) - t0

      # calculate weights
      t0 <- as.numeric(Sys.time())
      form <- c(
        "age"["age" %in% colnames(x)], "sex"["sex" %in% colnames(x)],
        variables
      ) |>
        paste0(collapse = " + ")
      form <- paste0("exposure ~ ", form)
      model <- glm(formula = form, data = x, family = binomial(link = "logit"))
      prob <- predict.glm(object = model, newdata = x, type = "response")
      x <- x |>
        mutate(ow = if_else(exposure == 1, 1 - prob, prob))
      time$weights <- time$weights + as.numeric(Sys.time()) - t0

      # fit outcome model
      t0 <- as.numeric(Sys.time())
      out <- glm(outcome ~ exposure, family = binomial(link = "logit"), data = x, weights = ow)
      log_out <- coef(out)[["exposure"]]
      se_out <- as_tibble(summary(out)$coefficients)[["Std. Error"]][2]
      result_nm$crude <- tibble(
        type = "crude",
        coef = log_out,
        coef_se = se_out
      )
      time$crude <- time$crude + as.numeric(Sys.time()) - t0

      # fit calibrate
      t0 <- as.numeric(Sys.time())
      log_or <- numeric()
      se_or <- numeric()
      for (k in seq_len(negative_control_outcomes)) {
        form <- sprintf("nco_%03i ~ exposure", k)
        out <- glm(formula = form, family = binomial(link = "logit"), data = x, weights = ow)
        sum <- as_tibble(summary(out)$coefficients)
        log_or[k] <- sum[["Estimate"]][2]
        se_or[k] <- sum[["Std. Error"]][2]
      }
      model <- fitSystematicErrorModel(
        logRr = log_or, seLogRr = se_or, trueLogRr = rep(0, length(log_or))
      )
      calibrated <- calibrateConfidenceInterval(
        logRr = log_out, seLogRr = se_out, model = model
      )
      result_nm$calibrated <- tibble(
        type = "calibrated",
        coef = calibrated$logRr,
        coef_se = calibrated$seLogRr
      )
      time$calibrated <- time$calibrated + as.numeric(Sys.time()) - t0

      # fit outcome with nco
      t0 <- as.numeric(Sys.time())
      form <- paste0(
        "outcome ~ exposure + ",
        paste0(sprintf("nco_%03i", seq_len(negative_control_outcomes)), collapse = " + ")
      )
      out <- glm(formula = form, family = binomial(link = "logit"), data = x, weights = ow)
      log_out <- coef(out)[["exposure"]]
      se_out <- as_tibble(summary(out)$coefficients)[["Std. Error"]][2]
      result_nm$adjusted <- tibble(
        type = "adjusted",
        coef = log_out,
        coef_se = se_out
      )
      time$adjusted <- time$adjusted + as.numeric(Sys.time()) - t0

      result_i[[nm]] <- bind_rows(result_nm) |>
        mutate(
          reference = "methotrexate",
          comparator = nm,
          coef_real = treatment_effect[[nm]],
          percentage_exposed = sum(x$exposed),
          unmeasured = length(unmeasured)
        )
    }

    # update result
    result[[i]] <- bind_rows(result_i) |>
      mutate(simulation_id = i)
  }

  # format result
  result <- result |>
    bind_rows() |>
    mutate(
      number_confounders = length(confounders_candidates),
      sample_size = sample_size,
      outcome_prev = outcome_prev,
      number_simulations = number_simulations,
      negative_control_outcomes = negative_control_outcomes
    )

  # time
  time$other <- as.numeric(Sys.time()) - ti - sum(unlist(time))
  cli_inform("{number_simulations} finished in: {round(sum(unlist(time)))} seconds:")
  imap(time, \(x, nm) cli_inform(c("*" = "{nm}: {round(x)} seconds."))) |>
    invisible()

  return(result)
}
fit_outcome_model <- function(data, number_confounders, unobserved_fraction) {
  cli_inform("Fitting initial model.")

  # select confounders
  confounders_candidates <- colnames(data) |>
    keep(\(x) startsWith(x = x, prefix = "var_")) |>
    sort()
  set.seed(123456)
  confounders_candidates <- c("age", "sex", sample(x = confounders_candidates, size = number_confounders - 2))
  outcome_formula <- paste0(
    "outcome ~ hydroxychloroquine + sulfasalazine + leflunomide + ",
    paste0(confounders_candidates, collapse = " + ")
  )

  # fit outcome model
  outcome_model <- glm(
    formula = outcome_formula,
    family = "binomial",
    data = data,
    control = glm.control(trace = TRUE)
  )

  list(outcome_model = outcome_model, confounders_candidates = confounders_candidates)
}
