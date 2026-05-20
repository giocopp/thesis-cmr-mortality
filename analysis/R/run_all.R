# ── Full pipeline ─────────────────────────────────────────────────────────────
# Reproduces every figure and table in output/ from raw data in data/raw/.
# See 01_download/README.md for raw-data acquisition.

renv::restore(prompt = FALSE)

options(width = 300, tibble.print_max = Inf, tibble.width = Inf)

root <- here::here("analysis", "R")

run_section <- function(folder) {
  files <- list.files(file.path(root, folder), pattern = "\\.R$", full.names = TRUE)
  purrr::walk(sort(files), source, echo = FALSE)
}

run_section("02_clean")
run_section("03_build")
run_section("04_descriptive")
run_section("05_analysis")
run_section("06_robustness")
