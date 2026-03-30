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
