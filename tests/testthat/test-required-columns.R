library(testthat)
required_cols <- c("sample_id", "assay", "project_name")

test_that("load_phenocycler_summary_files returns required columns", {
  testthat::local_mocked_bindings(
    .find_and_load_phenocycler_summary_files = function(data_dir = NULL) {
      data.frame(Sample = "A_1", `Image Tag` = "tag.tiff", stringsAsFactors = FALSE, check.names = FALSE)
    },
    load_meta_data = function() {
      data.frame(sample_id = "A_1", EBER_status = "positive", stringsAsFactors = FALSE)
    },
    .package = "EBVhelpR"
  )

  out <- load_phenocycler_summary_files()
  expect_s3_class(out, "data.frame")
  expect_true(all(required_cols %in% colnames(out)))
})

test_that("load_rnascope_summary_files returns required columns", {
  testthat::local_mocked_bindings(
    .find_and_load_rnascope_summary_files = function(data_dir = NULL) {
      data.frame(
        assay = EBV_ASSAY_TYPES$RNAScope_4plex,
        Sample = "A_1",
        SampleNumber = "A1",
        `Image Tag` = "tag.tiff", stringsAsFactors = FALSE, check.names = FALSE
      )
    },
    load_meta_data = function() {
      data.frame(sample_id = "A_1", EBER_status = "positive", stringsAsFactors = FALSE)
    },
    .package = "EBVhelpR"
  )

  out <- load_rnascope_summary_files()
  expect_s3_class(out, "data.frame")
  expect_true(all(required_cols %in% colnames(out)))
})

test_that("load_cell_source_files returns required columns", {
  tmp_dir <- withr::local_tempdir()
  dir.create(file.path(tmp_dir, "Phenocycler"), recursive = TRUE, showWarnings = FALSE)
  file.create(file.path(tmp_dir, "Phenocycler", "A_1.cell_data.csv"))
  file.create(file.path(tmp_dir, "Phenocycler", "CellPelletSlide_grp_PosCTL.cell_data.csv"))

  testthat::local_mocked_bindings(
    get_wrangled_cell_data_dir = function() tmp_dir,
    load_meta_data = function() {
      data.frame(sample_id = "A_1", EBER_status = "positive", stringsAsFactors = FALSE)
    },
    .package = "EBVhelpR"
  )

  out <- suppressWarnings(load_cell_source_files())
  expect_s3_class(out, "data.frame")
  expect_true(all(required_cols %in% colnames(out)))
})

test_that("get_tiff_file_path_df returns required columns", {
  out <- suppressWarnings(EBVhelpR:::get_tiff_file_path_df())
  expect_s3_class(out, "data.frame")
  expect_true(all(required_cols %in% colnames(out)))
})
