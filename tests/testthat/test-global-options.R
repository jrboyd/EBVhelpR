test_that("EBV_ASSAY_TYPES contains expected named assay labels", {
  expect_named(
    EBV_ASSAY_TYPES,
    c("phenocycler", "rnascope_4plex", "rnascope_3plex+IF")
  )
  expect_identical(EBV_ASSAY_TYPES$phenocycler, "Phenocycler")
  expect_identical(EBV_ASSAY_TYPES$rnascope_4plex, "RNAScope_4plex")
  expect_identical(EBV_ASSAY_TYPES[["rnascope_3plex+IF"]], "RNAScope_3plex+IF")
})
