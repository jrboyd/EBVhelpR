required_cols <- c("sample_id", "assay", "project_name")

expect_required_layout_live <- function(x) {
  expect_s3_class(x, "data.frame")
  expect_true(all(required_cols %in% colnames(x)))
  expect_type(x$sample_id, "character")
  expect_type(x$assay, "character")
  expect_type(x$project_name, "character")
}

has_original_data <- function() {
  tryCatch(dir.exists(get_original_cell_data_dir()), error = function(e) FALSE)
}

has_wrangled_cell_data <- function() {
  tryCatch({
    d <- get_wrangled_cell_data_dir()
    dir.exists(d) && length(list.files(d, pattern = "\\.cell_data\\.csv$", recursive = TRUE)) > 0
  }, error = function(e) FALSE)
}

test_that("[live] load_phenocycler_summary_files returns required layout", {
  skip_if_not(has_original_data(), "Original data directory not available on this machine.")

  data_dir <- get_original_cell_data_dir()
  has_pheno <- length(list.files(
    data_dir,
    pattern = "Summary.+csv",
    recursive = TRUE,
    full.names = TRUE
  )) > 0
  skip_if_not(has_pheno, "No phenocycler summary files found in data directory.")

  out <- load_phenocycler_summary_files()
  expect_required_layout_live(out)
})

test_that("[live] load_rnascope_summary_files returns required layout", {
  skip_if_not(has_original_data(), "Original data directory not available on this machine.")

  data_dir <- get_original_cell_data_dir()
  rscope_candidates <- list.files(
    data_dir,
    pattern = "RNA.+csv",
    recursive = TRUE,
    full.names = TRUE
  )
  rscope_candidates <- rscope_candidates[!grepl("Wide|Object", rscope_candidates)]
  rscope_candidates <- rscope_candidates[
    grepl("RNAScopeIF_Coexpression_2026-02-11", rscope_candidates) |
      grepl("RNAScope_Coexpression_2026-01-12.csv", rscope_candidates)
  ]
  skip_if_not(length(rscope_candidates) > 0, "No matching RNAscope summary files found.")

  out <- load_rnascope_summary_files()
  expect_required_layout_live(out)
})

test_that("[live] load_cell_source_files returns required layout", {
  skip_if_not(has_wrangled_cell_data(), "No wrangled .cell_data.csv files found.")

  out <- suppressWarnings(load_cell_source_files())
  expect_required_layout_live(out)
})

test_that("[live] get_tiff_file_path_df returns required layout", {
  out <- suppressWarnings(EBVhelpR:::get_tiff_file_path_df())
  expect_required_layout_live(out)
})
