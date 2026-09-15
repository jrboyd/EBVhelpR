library(testthat)

.gate_frame <- function(cd4_pos = 10, cd4_neg = 90, cd8_pos = 30, cd8_neg = 70) {
  df <- data.frame(
    `Image Tag` = "D-EB-10.qptiff",
    `CD4+TCells Cells` = cd4_pos,
    `CD4-Tcells Cells` = cd4_neg,
    `CD8+Tcells Cells` = cd8_pos,
    `CD8-Tcells Cells` = cd8_neg,
    `PD1+Tcells Cells` = 5,
    `PDL1+Bcells Cells` = 7,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  pct <- c(
    "% CD4+TCells Positive Cells", "% CD4-Tcells Positive Cells",
    "% CD8+TCells Positive Cells", "% CD8-Tcells Positive Cells",
    "% PD1+Tcells Positive Cells", "% PDL1+Bcells Positive Cells"
  )
  # names must match the renamer exactly, including the CD8 capitalisation
  pct[3] <- "% CD8+Tcells Positive Cells"
  for (p in pct) df[[p]] <- 1
  df
}

test_that("gate columns lose + and - on import", {
  out <- EBVhelpR:::.rename_phenocycler_gate_columns(.gate_frame())
  cn <- colnames(out)

  expect_false(any(grepl("[+-]", cn)))
  expect_true("CD4_pos_Tcells Cells" %in% cn)
  expect_true("CD4_neg_Tcells Cells" %in% cn)
  expect_true("CD8_pos_Tcells Cells" %in% cn)
  expect_true("CD8_neg_Tcells Cells" %in% cn)
  expect_true("PDL1_pos_Pax5_Bcells Cells" %in% cn)
  expect_true("% CD4_pos_Tcells Positive Cells" %in% cn)
})

test_that("renamed columns survive make.names() without colliding", {
  out <- EBVhelpR:::.rename_phenocycler_gate_columns(.gate_frame())
  mangled <- make.names(colnames(out), unique = FALSE)
  expect_equal(anyDuplicated(mangled), 0L)
})

test_that("a missing gate column is an error, not a silent pass", {
  df <- .gate_frame()
  df[["CD4+TCells Cells"]] <- NULL
  expect_error(
    EBVhelpR:::.rename_phenocycler_gate_columns(df),
    "missing expected gate column"
  )
})

test_that("the CD4/CD8 complement assertion passes on balanced data", {
  out <- EBVhelpR:::.rename_phenocycler_gate_columns(.gate_frame())
  expect_silent(EBVhelpR:::.assert_phenocycler_gate_complement(out))
})

test_that("the CD4/CD8 complement assertion fires when the gates disagree", {
  out <- EBVhelpR:::.rename_phenocycler_gate_columns(
    .gate_frame(cd4_pos = 10, cd4_neg = 90, cd8_pos = 30, cd8_neg = 71)
  )
  expect_error(
    EBVhelpR:::.assert_phenocycler_gate_complement(out),
    "do not partition the same"
  )
})
