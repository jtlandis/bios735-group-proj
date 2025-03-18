data_url <- url(
  "https://data.mendeley.com/public-files/datasets/njdkntcpc9/files/08bb4f43-6dfa-4995-b268-42fa0690ba6b/file_downloaded"
)

data_set_raw <- readr::read_csv(data_url)

data_set_tidy <- tidyr::pivot_longer(
  data_set_raw,
  cols = -DATE,
  names_to = c("type", "brand", "item"),
  names_pattern = "([^_]+)_([^_]+)_([^_]+)",
  values_to = "count"
) |>
  tidyr::pivot_wider(
    names_from = type,
    values_from = count
  )

usethis::use_data(data_set_raw, overwrite = TRUE)
usethis::use_data(data_set_tidy, overwrite = TRUE)
