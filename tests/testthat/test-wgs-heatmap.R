library(testthat)

.pileup <- function(samples = c("S1", "S2"), status = "Positive") {
  expand.grid(sample = samples, x = c(0, 100), stringsAsFactors = FALSE) |>
    transform(y = 1, EBER_status = status)
}

.count_summary <- function(levels_in) {
  data.frame(sample_id = factor(levels_in, levels = levels_in))
}

test_that("samples missing from wgs_count_summary warn and are dropped", {
  pileup <- .pileup(c("S1", "S2", "S3"))
  expect_warning(
    p <- plot_wgs_pileup_heatmap(pileup, .count_summary(c("S1", "S2"))),
    "absent from `wgs_count_summary`"
  )
  # The dropped sample must not survive as an NA row collecting every
  # out-of-level sample at one y position.
  expect_false(any(is.na(p$data$sample)))
  expect_setequal(as.character(unique(p$data$sample)), c("S1", "S2"))
})

test_that("no warning when every sample is in the level set", {
  pileup <- .pileup(c("S1", "S2"))
  expect_no_warning(plot_wgs_pileup_heatmap(pileup, .count_summary(c("S1", "S2"))))
})

test_that("an EBER_status with no color is an error, not a transparent bar", {
  pileup <- .pileup(c("S1", "S2"), status = "need info")
  expect_error(
    plot_wgs_pileup_heatmap(
      pileup,
      .count_summary(c("S1", "S2")),
      status_colors = c(Positive = "red", Negative = "blue")
    ),
    "No color supplied for EBER_status"
  )
})

test_that("the status bar still renders without wgs_count_summary", {
  # `sample` is character here, so as.numeric() used to yield all NA and the
  # annotation silently vanished.
  p <- plot_wgs_pileup_heatmap(.pileup(c("S1", "S2")))
  expect_s3_class(p, "ggplot")
  expect_length(p$layers, 2L)
})
