library(testthat)

mock_rnascope <- function(sample_numbers) {
  testthat::local_mocked_bindings(
    .find_and_load_rnascope_summary_files = function(data_dir = NULL) {
      data.frame(
        assay = EBV_ASSAY_TYPES$RNAScope_4plex,
        Sample = NA_character_,
        SampleNumber = sample_numbers,
        combo = "EBER1", count = 10,
        stringsAsFactors = FALSE
      )
    },
    load_meta_data = function() {
      data.frame(sample_id = "CTEBV_15", EBER_status = "Positive", stringsAsFactors = FALSE)
    },
    .package = "EBVhelpR",
    .env = parent.frame()
  )
}

test_that("load_rnascope_summary_files warns when probe controls share a sample_id", {
  mock_rnascope(c("CTEBV15", "CTEBV15PosCTL", "CTEBV15NegCTL"))
  expect_warning(out <- load_rnascope_summary_files(), "probe-control sections")
  expect_equal(out$sample_id, rep("CTEBV_15", 3))
  expect_setequal(out$probe_control, c("", "positive_probe", "negative_probe"))
})

test_that("load_rnascope_summary_files is silent without probe controls", {
  mock_rnascope(c("CTEBV15"))
  expect_no_warning(load_rnascope_summary_files())
})

test_that("the probe-control warning can be switched off", {
  mock_rnascope(c("CTEBV15", "CTEBV15PosCTL"))
  withr::local_options(EBVhelpR.warn_probe_controls = FALSE)
  expect_no_warning(load_rnascope_summary_files())
})

test_that("only the PosCTL / NegCTL suffixes mark a probe control", {
  mock_rnascope(c("Posterior1", "Negative2", "CTEBV15NegCTL"))
  out <- suppressWarnings(load_rnascope_summary_files())
  expect_equal(out$probe_control[out$sample_id == "Posterior1"], "")
  expect_equal(out$probe_control[out$sample_id == "Negative2"], "")
  expect_equal(out$probe_control[out$sample_id == "CTEBV_15"], "negative_probe")
})

test_that(".warn_probe_controls ignores tables without a probe_control column", {
  expect_no_warning(EBVhelpR:::.warn_probe_controls(data.frame(sample_id = "A"), "f"))
})

test_that("check_summary_against_cell_files flags a summary that pools extra sections", {
  a <- EBV_ASSAY_TYPES$RNAScope_4plex
  summary_df <- data.frame(
    assay = a,
    sample_id = c("GM12878", "GM12878", "CTEBV_15", "CTEBV_15"),
    probe_control = c("", "", "", "positive_probe"),
    combo = c("EBER1", "Unstained_Cells", "EBER1", "EBER1"),
    count = c(60, 90, 40, 70),              # GM12878 summary holds 150 cells
    stringsAsFactors = FALSE
  )
  cell_files_df <- data.frame(
    assay = a,
    sample_id = c("GM12878", "CTEBV_15", "CTEBV_15"),
    probe_control = c("", "", "positive_probe"),
    cell_count = c(101, 41, 71),            # wc -l: one header line each
    stringsAsFactors = FALSE
  )
  expect_warning(chk <- check_summary_against_cell_files(summary_df, cell_files_df),
                 "GM12878")
  expect_false(chk$match[chk$unique_id == "GM12878"])
  expect_equal(chk$n_summary[chk$unique_id == "GM12878"], 150)
  expect_equal(chk$n_cell_files[chk$unique_id == "GM12878"], 100)
  expect_true(chk$match[chk$unique_id == "CTEBV_15"])
  expect_true(chk$match[chk$unique_id == "CTEBV_15 positive_probe"])

  ok <- summary_df[summary_df$sample_id == "CTEBV_15", ]
  expect_no_warning(check_summary_against_cell_files(ok, cell_files_df[-1, ]))
})

test_that("check_summary_against_cell_files reports missing columns", {
  expect_error(check_summary_against_cell_files(data.frame(assay = "x"), data.frame()),
               "summary_df lacks")
})
