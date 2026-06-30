# read data
data <- readr::read_csv(file = here::here("clean_data.csv"), show_col_types = FALSE)

covariates <- colnames(data) |>
  purrr::keep(\(x) startsWith(x = x, prefix = "var_"))

# check exposure balance
ref <- "methotrexate"
comp <- data$exposure |>
  unique() |>
  purrr::keep(\(x) x != ref)
numeric <- "var_0001"
binary <- covariates[covariates != numeric]

# numeric data
est <- data |>
  dplyr::select(dplyr::all_of(c("exposure", numeric))) |>
  tidyr::pivot_longer(cols = dplyr::all_of(numeric), names_to = "variable", values_to = "value") |>
  dplyr::group_by(exposure, variable) |>
  dplyr::summarise(mean = mean(value), sd = sd(value), .groups = "drop")
smds_numeric <- est |>
  dplyr::filter(exposure == ref) |>
  dplyr::rename(reference = "exposure", mean_ref = "mean", sd_ref = "sd") |>
  dplyr::inner_join(
    est |>
      dplyr::filter(exposure %in% comp) |>
      dplyr::rename(comparator = "exposure", mean_comp = "mean", sd_comp = "sd"),
    by = "variable",
    relationship = "one-to-many"
  ) |>
  dplyr::mutate(
    value_ref = sprintf("%.1f (%.1f)", mean_ref, sd_ref),
    value_comp = sprintf("%.1f (%.1f)", mean_comp, sd_comp),
    smd = (mean_comp - mean_ref) / sqrt(sd_comp ** 2 + sd_ref ** 2)
  ) |>
  dplyr::select(reference, comparator, variable, value_ref, value_comp, smd)

# binary variables
est <- data |>
  dplyr::select(dplyr::all_of(c("exposure", binary))) |>
  dplyr::group_by(exposure) |>
  dplyr::summarise(dplyr::across(.cols = dplyr::all_of(binary), .fns = list(p = \(x) mean(x)), .names = "{.col}")) |>
  tidyr::pivot_longer(cols = dplyr::all_of(binary), names_to = "variable", values_to = "p")
smds_binary <- est |>
  dplyr::filter(exposure == ref) |>
  dplyr::rename(reference = "exposure", p_ref = "p") |>
  dplyr::inner_join(
    est |>
      dplyr::filter(exposure %in% comp) |>
      dplyr::rename(comparator = "exposure", p_comp = "p"),
    by = "variable",
    relationship = "one-to-many"
  ) |>
  dplyr::mutate(
    smd = (p_comp - p_ref) / sqrt((p_comp * (1 - p_comp) + p_ref * (1 - p_ref)) / 2),
    value_ref = sprintf("%.1f%%", 100 * p_ref),
    value_comp = sprintf("%.1f%%", 100 * p_comp)
  ) |>
  dplyr::select(reference, comparator, variable, value_ref, value_comp, smd)

# smds
smds <- dplyr::union_all(smds_numeric, smds_binary)

# visualise
x <- smds |>
  dplyr::group_by(reference, comparator) |>
  dplyr::group_split() |>
  purrr::map(\(x) {
    x |>
      dplyr::filter(smd > 0) |>
      dplyr::arrange(dplyr::desc(smd)) |>
      dplyr::mutate(x = dplyr::row_number() - 0.5) |>
      dplyr::union_all(
        x |>
          dplyr::filter(smd < 0) |>
          dplyr::arrange(smd) |>
          dplyr::mutate(x = - dplyr::row_number() + 0.5, smd = abs(smd))
      )
  }) |>
  dplyr::bind_rows()
ggplot2::ggplot(data = x) +
  ggplot2::geom_point(mapping = ggplot2::aes(x = x, y = smd, colour = comparator)) +
  ggplot2::geom_hline(yintercept = 0.1) +
  ggplot2::coord_cartesian(xlim = c(-100, 100), ylim = c(0, 3))

# concept frequency
data |>
  dplyr::summarise(dplyr::across(
    dplyr::starts_with("var_"),
    mean,
    .names = "{.col}"
  )) |>
  tidyr::pivot_longer(dplyr::everything()) |>
  dplyr::filter(name != "var_0001") |>
  dplyr::arrange(dplyr::desc(value)) |>
  dplyr::mutate(x = dplyr::row_number()) |>
  ggplot2::ggplot(mapping = ggplot2::aes(x = x, y = value)) +
  ggplot2::geom_point()
