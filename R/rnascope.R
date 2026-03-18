#' Load RNAscope summary CSV files
#'
#' Reads RNAscope coexpression summary files from the configured data directory
#' and combines selected assay versions into one data frame.
#'
#' @param data_dir Optional directory containing source summary files. If `NULL`,
#'   the package resolver checks `EBVHELPER_DATA_DIR` and known default paths.
#'
#' @return A data frame with one row per record and a `source` column indicating
#'   the assay source.
#' @export
load_rnascope_summary_files <- function(data_dir = NULL) {
  if (is.null(data_dir)) {
    data_dir <- .get_data_dir()
  }
  stopifnot(dir.exists(data_dir))

  res_files <- list.files(
    data_dir,
    pattern = "RNA.+csv",
    recursive = TRUE,
    full.names = TRUE
  )

  res_files <- res_files[!grepl("Wide", res_files)]
  res_files <- res_files[!grepl("Object", res_files)]
  res_files <- res_files[
    grepl("RNAScopeIF_Coexpression_2026-02-11", res_files) |
      grepl("RNAScope_Coexpression_2026-01-12.csv", res_files)
  ]

  if (!length(res_files)) {
    stop("No RNAscope summary files found for configured patterns.", call. = FALSE)
  }

  file_rename <- c(
    "RNAScope_Coexpression_2026-01-12.csv" = "4plex",
    "RNAScopeIF_Coexpression_2026-02-11.csv" = "3plex+IF"
  )

  names(res_files) <- file_rename[basename(res_files)]
  keep <- !is.na(names(res_files))
  res_files <- as.list(res_files[keep])

  if (!length(res_files)) {
    stop("No RNAscope files matched known filenames.", call. = FALSE)
  }

  all_dt_l <- .load_csv_list(res_files)
  dplyr::bind_rows(all_dt_l, .id = "source")
}

#' Harmonize RNAscope summaries with EBV status metadata
#'
#' Loads RNAscope summary data and metadata, maps sample_id identifiers, and appends
#' `EBER_status` annotation.
#'
#' @param data_dir Optional directory containing source summary files. If `NULL`,
#'   the package resolver checks `EBVHELPER_DATA_DIR` and known default paths.
#'
#' @return A data frame with harmonized identifiers and EBV status annotation.
#' @export
harmonize_rnascope_summary_files <- function(data_dir = NULL) {
  rscope_dt <- load_rnascope_summary_files(data_dir = data_dir)
  meta_df <- load_meta_data()

  s2s <- meta_df$sample_id
  names(s2s) <- meta_df$SampleStripped
  rscope_dt$SampleID <- s2s[rscope_dt$SampleNumber]

  anno_df <- dplyr::select(meta_df, SampleID = .data$sample_id, .data$EBER_status)
  rscope_dt <- merge(rscope_dt, anno_df, all.x = TRUE)
  rscope_dt <- dplyr::mutate(
    rscope_dt,
    EBER_status = ifelse(is.na(.data$EBER_status), "need info", .data$EBER_status)
  )

  rscope_dt
}