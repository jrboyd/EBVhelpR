#' Load phenocycler summary CSV files
#'
#' Reads all phenocycler summary CSV files from the configured data directory,
#' combines them, and derives a standardized `Sample` identifier.
#'
#' @param data_dir Optional directory containing source summary files. If `NULL`,
#'   the package resolver checks `EBVHELPER_DATA_DIR` and known default paths.
#'
#' @return A data frame with one row per record and a `source` column indicating
#'   the file group.
#' @export
load_phenocycler_summary_files <- function(data_dir = NULL) {
  if (is.null(data_dir)) {
    data_dir <- .get_data_dir()
  }
  stopifnot(dir.exists(data_dir))

  res_files <- list.files(
    data_dir,
    pattern = "Summary.+csv",
    recursive = TRUE,
    full.names = TRUE
  )

  if (!length(res_files)) {
    stop("No phenocycler summary files found.", call. = FALSE)
  }

  names(res_files) <- basename(dirname(res_files))
  all_dt_l <- .load_csv_list(as.list(res_files))

  dt <- dplyr::bind_rows(all_dt_l, .id = "source")
  if (!"Image Tag" %in% colnames(dt)) {
    stop("Expected column `Image Tag` was not found in phenocycler summaries.", call. = FALSE)
  }

  dt <- dt |>
    dplyr::mutate(Sample = sub("\\..+", "", .data$`Image Tag`)) |>
    dplyr::mutate(Sample = sub("_Scan.+", "", .data$Sample)) |>
    dplyr::mutate(Sample = gsub("-", "", .data$Sample))

  dt
}

#' Harmonize phenocycler summaries with EBV status metadata
#'
#' Loads phenocycler summary data and metadata, maps sample identifiers, and
#' appends `EBV_Status` annotation.
#'
#' @param data_dir Optional directory containing source summary files. If `NULL`,
#'   the package resolver checks `EBVHELPER_DATA_DIR` and known default paths.
#'
#' @return A data frame with harmonized identifiers and EBV status annotation.
#' @export
harmonize_phenocycler_summary_files <- function(data_dir = NULL) {
  pcycler_dt <- load_phenocycler_summary_files(data_dir = data_dir)
  meta_df <- load_meta_data(data_dir = data_dir)

  s2s <- meta_df$sample
  names(s2s) <- meta_df$SampleStripped
  pcycler_dt$SampleID <- s2s[pcycler_dt$Sample]

  anno_df <- dplyr::select(meta_df, SampleID = .data$sample, .data$EBV_Status)
  pcycler_dt <- merge(pcycler_dt, anno_df, all.x = TRUE)
  pcycler_dt <- dplyr::mutate(
    pcycler_dt,
    EBV_Status = ifelse(is.na(.data$EBV_Status), "need info", .data$EBV_Status)
  )

  pcycler_dt
}