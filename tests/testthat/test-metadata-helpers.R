test_that("expand helpers split composite sample IDs", {
  df <- data.frame(
    sample_id = c("A_1", "B_2 previously C_3", "D_4/E_5"),
    EBER_status = c("positive", "negative", "unknown"),
    stringsAsFactors = FALSE
  )

  previous_out <- EBVhelpR:::.expand_previous_samples(df)
  expect_true(all(c("B_2", "C_3") %in% previous_out$sample_id))

  slash_out <- EBVhelpR:::.expand_slash_samples(df)
  expect_true(all(c("D_4", "E_5") %in% slash_out$sample_id))
})

test_that("load_meta_data rejects a sheet too narrow to name three columns", {
  # The check used to require only two columns while assigning three names, so a
  # narrow sheet failed one line later with an opaque 'names' attribute error.
  narrow <- file.path(withr::local_tempdir(), "narrow.xlsx")
  openxlsx::write.xlsx(
    data.frame(Sample = "D-EB-11", Status = "Positive", stringsAsFactors = FALSE),
    narrow
  )

  testthat::local_mocked_bindings(
    .get_status_file = function() narrow,
    .package = "EBVhelpR"
  )

  expect_error(load_meta_data(), "at least three columns")
})

test_that("load_meta_data drops Notes before counting columns", {
  wide <- file.path(withr::local_tempdir(), "wide.xlsx")
  openxlsx::write.xlsx(
    data.frame(
      Sample = "D-EB-11",
      Status = "Positive",
      Type = "biopsy",
      Notes = "reclassified 2026-09-15",
      stringsAsFactors = FALSE
    ),
    wide
  )

  testthat::local_mocked_bindings(
    .get_status_file = function() wide,
    .package = "EBVhelpR"
  )

  out <- load_meta_data()
  expect_false("Notes" %in% colnames(out))
  expect_equal(colnames(out)[1:3], c("sample_id", "EBER_status", "sample_type"))
  expect_equal(out$sample_id, "D_EB_11")
})
