# ── Full pipeline ────────────────────────────────────────────────────────────

renv::restore(prompt = FALSE)

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
