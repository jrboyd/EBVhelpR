library(testthat)

# Build a CellQueryInfo directly, bypassing CellQuery()'s data loading.
# `D_EB_33` is present in summary_df but absent from all_cell_files_df, which is
# the real situation created by dropping the multi-sample D_EB_12_D_EB_33 image.
.make_query <- function() {
  summary_df <- EBVhelpR:::.df_prep(data.frame(
    sample_id = c("D_EB_10", "D_EB_33", "CTEBV_1"),
    probe_control = c("", "", "negative_probe"),
    stringsAsFactors = FALSE
  ))
  cell_files_df <- EBVhelpR:::.df_prep(data.frame(
    sample_id = c("D_EB_10", "CTEBV_1"),
    probe_control = c("", "negative_probe"),
    stringsAsFactors = FALSE
  ))
  tiff_df <- EBVhelpR:::.df_prep(data.frame(
    sample_id = c("D_EB_10", "D_EB_33"),
    probe_control = c("", ""),
    stringsAsFactors = FALSE
  ))

  methods::new(
    "CellQueryInfo",
    summary_df = summary_df,
    all_cell_files_df = cell_files_df,
    meta_data_df = data.frame(
      sample_id = c("D_EB_10", "D_EB_33", "CTEBV_1"),
      EBER_status = "Positive",
      stringsAsFactors = FALSE
    ),
    selected_sample_ids = unique(cell_files_df$sample_id),
    selected_unique_ids = unique(cell_files_df$unique_id),
    tiff_paths_df = tiff_df,
    assay_type = "Phenocycler"
  )
}

test_that(".df_prep builds unique_id from sample_id and probe_control", {
  q <- .make_query()
  expect_equal(
    q@summary_df$unique_id,
    c("D_EB_10", "D_EB_33", "CTEBV_1 negative_probe")
  )
})

test_that("set_selected_unique_ids derives selected_sample_ids", {
  q <- .make_query()
  q <- set_selected_unique_ids(q, c("D_EB_10", "CTEBV_1 negative_probe"))

  expect_setequal(q@selected_unique_ids, c("D_EB_10", "CTEBV_1 negative_probe"))
  expect_setequal(q@selected_sample_ids, c("D_EB_10", "CTEBV_1"))
})

test_that("set_selected_sample_ids derives selected_unique_ids", {
  q <- .make_query()
  q <- set_selected_sample_ids(q, "CTEBV_1")

  expect_setequal(q@selected_unique_ids, "CTEBV_1 negative_probe")
  expect_setequal(q@selected_sample_ids, "CTEBV_1")
})

test_that("a summary-only sample survives selection (issue 1: n=29 vs n=30)", {
  q <- .make_query()
  # Select everything the summary knows about, as an overlap/upset step would.
  q <- set_selected_unique_ids(q, q@summary_df$unique_id)

  expect_true("D_EB_33" %in% q@selected_unique_ids)
  expect_true("D_EB_33" %in% get_query_summary_df(q)$sample_id)
  expect_equal(nrow(get_query_summary_df(q)), 3L)
})

test_that("probe-control ids are selectable (issue 2)", {
  q <- .make_query()
  # Before the fix selected_unique_ids held bare sample ids, so a probe-control
  # unique_id could never be selected.
  expect_true("CTEBV_1 negative_probe" %in% q@selected_unique_ids)

  q <- set_selected_unique_ids(q, "CTEBV_1 negative_probe")
  expect_equal(nrow(get_query_summary_df(q)), 1L)
  expect_equal(get_query_summary_df(q)$sample_id, "CTEBV_1")
})

test_that("unknown ids warn and are dropped", {
  q <- .make_query()
  expect_warning(
    q2 <- set_selected_unique_ids(q, c("D_EB_10", "NOT_A_SAMPLE")),
    "not present in this query"
  )
  expect_setequal(q2@selected_unique_ids, "D_EB_10")

  expect_warning(
    q3 <- set_selected_sample_ids(q, c("D_EB_10", "NOT_A_SAMPLE")),
    "not present in this query"
  )
  expect_setequal(q3@selected_sample_ids, "D_EB_10")
})

test_that("all three getters agree on a selection", {
  q <- .make_query()
  q <- set_selected_unique_ids(q, "D_EB_10")

  expect_equal(unique(get_query_summary_df(q)$unique_id), "D_EB_10")
  expect_equal(unique(get_query_cell_files_df(q)$unique_id), "D_EB_10")
  expect_equal(unique(get_query_tiff_paths_df(q)$unique_id), "D_EB_10")
})

test_that("get_query_cell_files_df honours selected_sample_ids", {
  q <- .make_query()
  q <- set_selected_unique_ids(q, "CTEBV_1 negative_probe")
  # Previously this getter applied the same unique_id filter twice and never
  # consulted selected_sample_ids, so it disagreed with the other two.
  expect_equal(nrow(get_query_cell_files_df(q)), 1L)
  expect_equal(get_query_cell_files_df(q)$sample_id, "CTEBV_1")
})
