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
.find_and_load_rnascope_summary_files <- function(data_dir = NULL) {
  if (is.null(data_dir)) {
    data_dir <- get_original_cell_data_dir()
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
  # res_files <- res_files[
  #   grepl("RNAScopeIF_Coexpression_2026-02-11", res_files) |
  #     grepl("RNAScope_Coexpression_2026-01-12.csv", res_files)
  # ]
  res_files <- res_files[
      !grepl("Sara", res_files)
  ]

  if (!length(res_files)) {
    stop("No RNAscope summary files found for configured patterns.", call. = FALSE)
  }

  file_rename <- c(
    "RNAScope_Coexpression_2026-01-12.csv" = EBV_ASSAY_TYPES$rnascope_4plex,
    "RNAScopeIF_Coexpression_2026-02-11.csv" = EBV_ASSAY_TYPES$`rnascope_3plex+IF`,
    "4PlexRNAScopeCellPelletSingleExpression_2026-01-12.csv" = EBV_ASSAY_TYPES$rnascope_4plex,
    "RNAScopeIFCellPelletSingleExpression_2026-02-11.csv" = EBV_ASSAY_TYPES$`rnascope_3plex+IF`
  )
  names(res_files) <- file_rename[basename(res_files)]
  #
  # keep <- !is.na(names(res_files))
  # res_files <- as.list(res_files[keep])

  if (!length(res_files)) {
    stop("No RNAscope files matched known filenames.", call. = FALSE)
  }

  all_dt_l <- .load_csv_list(res_files)

  dplyr::bind_rows(all_dt_l, .id = "assay")
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
#' @examples
#' \dontrun{
#' rscope_df <- load_rnascope_summary_files()
#' }
#' @export
load_rnascope_summary_files <- function(data_dir = NULL) {
  rscope_dt <- .find_and_load_rnascope_summary_files(data_dir = data_dir)
  meta_df <- load_meta_data()

  # rscope_dt = rscope_dt %>% dplyr::mutate(SampleStripped := sub("(NegCTL)|(PosCTL)", "", SampleNumber))
  rscope_dt = rscope_dt %>%
      dplyr::mutate(sample_id = ifelse(is.na(SampleNumber), Sample, SampleNumber)) %>%
      dplyr::mutate(sample_id = ifelse(is.na(SampleNumber), sample_id, sub("CTEBV", "CTEBV_", sample_id))) %>%
      dplyr::mutate(sample_id = ifelse(is.na(SampleNumber), sample_id, sub("DEB", "D_EB_", sample_id))) %>%
      dplyr::mutate(sample_id = sub("GM185022$", "GM18502", sample_id))


  rscope_dt$probe_control = "N/A"
  rscope_dt = rscope_dt %>% dplyr::mutate(probe_control = ifelse(grepl("[Nn]eg", Sample) | grepl("[Nn]eg", SampleNumber), "negative_probe", probe_control))
  rscope_dt = rscope_dt %>% dplyr::mutate(probe_control = ifelse(grepl("[Pp]os", Sample) | grepl("[Pp]os", SampleNumber), "positive_probe", probe_control))
  rscope_dt = rscope_dt %>% dplyr::mutate(sample_id = sub("NegCTL", "", sample_id))  %>% dplyr::mutate(sample_id = sub("PosCTL", "", sample_id))

  setdiff(rscope_dt$sample_id, meta_df$sample_id)
  setdiff(meta_df$sample_id, rscope_dt$sample_id)

  rscope_dt <- merge(rscope_dt, meta_df, all.x = TRUE, by = "sample_id")
  rscope_dt <- dplyr::mutate(
    rscope_dt,
    EBER_status = ifelse(is.na(.data$EBER_status), "need info", .data$EBER_status)
  )

  if (!"project_name" %in% colnames(rscope_dt)) {
    rscope_dt$project_name <- assay_to_project_name[rscope_dt$assay]
  }

  rscope_dt
}


#' Load cell-level source data files
#'
#' Discovers per-sample cell data CSVs written by
#' [write_all_package_data()] and assembles them into a
#' single data frame. Cohort samples are merged with clinical metadata from
#' [load_meta_data()]; control samples have their identifiers
#' standardized.
#'
#' @return A data frame with all cell source file paths, sample identifiers,
#'   assay grouping, and (for cohort samples) clinical metadata.
#' @examples
#' \dontrun{
#' cell_df <- load_cell_source_files()
#' }
#' @export
load_cell_source_files <- function() {
  pkg_data_dir <- get_wrangled_cell_data_dir()

  all_cell_data_files <- list.files(
    pkg_data_dir,
    pattern = "\\.cell_data\\.csv$",
    recursive = TRUE
  )

  cell_df <- data.frame(file = all_cell_data_files, stringsAsFactors = FALSE)
  cell_df <- tidyr::separate(
    cell_df,
    col = "file",
    sep = "/",
    into = c("assay", "image_name"),
    remove = FALSE
  )
  cell_df$image_name <- sub("\\..+", "", cell_df$image_name)
  cell_df$name <- cell_df$image_name

  cell_df <- cell_df |>
    dplyr::mutate(
      name = sub("_ ?[cC]ro.+", "", .data$name),
      name = sub("^[0-9]+_", "", .data$name),
      name = sub("EBER-LMP1-EBNA1", "PosCTL", .data$name),
      name = gsub("_Rescanned", "", .data$name)
    )  |>
    dplyr::mutate( #phenocycler specific
      name = sub("_Control_", "_NegCTL_", .data$name),
      name = sub("_Test_", "_PosCTL_", .data$name)
    ) |>
    dplyr::mutate(
      sample_type = "cohort",
      sample_type = ifelse(grepl("CTL", .data$name), "control", .data$sample_type)
    ) |>
    dplyr::mutate(
      name = gsub("-", "_", .data$name),
      name = sub("^CTEBV", "CTEBV_", .data$name),
      name = sub("_Scan[0-9]", "", .data$name),
      name = sub("CellPellet_", "", .data$name)
    )

  cell_df$sample_id <- cell_df$name

  cell_df.by_type <- split(cell_df, cell_df$sample_type)
  meta_df <- load_meta_data()

  cell_df.by_type$cohort <- merge(cell_df.by_type$cohort, meta_df, all.x = TRUE)

  warning("removing D_EB_12_D_EB_33_Scan5, multiple samples in same image")
  cell_df.by_type$cohort = cell_df.by_type$cohort %>% subset(!grepl("D_EB_12_D_EB_33", name))

  valid = cell_df.by_type$cohort$sample_id %in% meta_df$sample_id
  cell_df.by_type$cohort[!valid,]

  stopifnot(all(cell_df.by_type$cohort$sample_id %in% meta_df$sample_id))

  cell_df.by_type$control <- cell_df.by_type$control |>
    dplyr::mutate(sample_id = sub("_2$", "", .data$sample_id)) |>
    dplyr::mutate(
      control_type = ifelse(grepl("PosCTL", .data$sample_id), "PosCTL", "NegCTL")
    ) |>
    dplyr::mutate(sample_id = sub("CellPelletSlide_", "", .data$sample_id))|>
    dplyr::mutate(sample_id = sub("_Phenocycler", "", .data$sample_id))
  cell_df.by_type$control = cell_df.by_type$control |>
    dplyr::group_by(.data$file) |>
    dplyr::mutate(
      control_group = sub(.data$control_type, "", .data$sample_id),
      control_group = sub("^_", "", .data$control_group),
      control_group = sub("_$", "", .data$control_group),
      sample_id = paste(.data$control_group, .data$control_type, sep = "_")
    ) |>
    dplyr::ungroup()

  final_df = dplyr::bind_rows(cell_df.by_type$cohort, cell_df.by_type$control)

  final_df$probe_control = "N/A"
  final_df = final_df %>% dplyr::mutate(probe_control = ifelse(grepl("[Nn]eg", sample_id), "negative_probe", probe_control))
  final_df = final_df %>% dplyr::mutate(probe_control = ifelse(grepl("[Pp]os", sample_id), "positive_probe", probe_control))

  final_df$sample_id <- sub("_?NegCTL", "", final_df$sample_id)
  final_df$sample_id <- sub("_?PosCTL", "", final_df$sample_id)
  final_df$name = NULL
  final_df$project_name <- assay_to_project_name[final_df$assay]
  final_df
}
