
## load raw data
dir <- tempdir()
zip::unzip(zipfile = here::here("data_sent.zip"), exdir = dir)
file <- file.path(dir, "data_sent.csv")
data <- readr::read_csv(file = file) |>
  dplyr::mutate(
    exposure = exposure |>
      stringr::str_split(pattern = "_") |>
      purrr::map_chr(\(x) x[2]),
    id = dplyr::row_number()
  ) |>
  dplyr::select(!"...1") |>
  dplyr::relocate("id", "exposure", "outcome", "age" = "var2", "sex" = "var3")
unlink(x = file)

# rename
cn <- colnames(data)
id <- 6:length(cn)
cn[id] <- sprintf("var_%04i", id - 5)
colnames(data) <- cn

# trim frequency
minFrequency <- 0.005
selectedCovariates <- data |>
  dplyr::summarise(dplyr::across(
    dplyr::starts_with("var_"),
    list(p = \(x) mean(x)),
    .names = "{.col}"
  )) |>
  tidyr::pivot_longer(dplyr::everything()) |>
  dplyr::filter(value >= minFrequency & 1 - value >= minFrequency) |>
  dplyr::pull("name")
data <- data |>
  dplyr::select(dplyr::all_of(c("id", "exposure", "outcome", "age", "sex", sort(selectedCovariates))))

# save clean data
readr::write_csv(x = data, file = here::here("clean_data.csv"))
