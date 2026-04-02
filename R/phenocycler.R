#' Load Phenocycler summary CSV files
#'
#' Reads all Phenocycler summary CSV files from the configured data directory,
#' combines them, and derives a standardized `Sample` identifier.
#'
#' @param data_dir Optional directory containing source summary files. If `NULL`,
#'   the package resolver checks `EBVHELPER_DATA_DIR` and known default paths.
#'
#' @return A data frame with one row per record and a `source` column indicating
#'   the file group.
.find_and_load_phenocycler_summary_files <- function(data_dir = NULL) {
    if (is.null(data_dir)) {
        data_dir <- get_original_cell_data_dir()
    }
    stopifnot(dir.exists(data_dir))

    res_files <- list.files(
        data_dir,
        pattern = "Summary.+csv",
        recursive = TRUE,
        full.names = TRUE
    )

    if (!length(res_files)) {
        stop("No Phenocycler summary files found.", call. = FALSE)
    }

    names(res_files) <- basename(dirname(res_files))
    all_dt_l <- .load_csv_list(as.list(res_files))

    dt <- dplyr::bind_rows(all_dt_l, .id = "source")
    if (!"Image Tag" %in% colnames(dt)) {
        stop("Expected column `Image Tag` was not found in Phenocycler summaries.", call. = FALSE)
    }

    # dt.cell_pellet = dt %>% dplyr::filter(grepl("CellPelletSlide", `Image Tag`))
    #
    # str_last = function(x){
    #     sapply(strsplit(sub("\\..+", "", x), "_"), function(xx){xx[length(xx)]})
    # }
    #
    # dt.cell_pellet = dt.cell_pellet %>% dplyr::mutate(Sample = paste(sep = "_",
    #                                                                  str_last(`Image Tag`),
    #                                                                  ifelse(grepl("control", "Image Tag"), "NegCTL", "PosCTL")
    # ))
    # dt.cell_pellet$Sample
    #
    # dt.main = dt %>% dplyr::filter(!grepl("CellPelletSlide", `Image Tag`))
    #
    # dt.main <- dt.main |>
    #     dplyr::mutate(Sample = sub("\\..+", "", .data$`Image Tag`)) |>
    #     dplyr::mutate(Sample = sub("_Scan.+", "", .data$Sample)) |>
    #     dplyr::mutate(Sample = gsub("-", "", .data$Sample))

    # dt = rbind(dt.main, dt.cell_pellet)
    dt
}

#' Harmonize Phenocycler summaries with EBV status metadata and ensure compatible sample ids.
#'
#' Loads Phenocycler summary data and metadata, maps sample identifiers, and
#' appends `EBER_status` annotation.
#'
#' @param data_dir Optional directory containing source summary files. If `NULL`,
#'   the package resolver checks `EBVHELPER_DATA_DIR` and known default paths.
#'
#' @return A data frame with harmonized identifiers and EBV status annotation.
#' @examples
#' \dontrun{
#' pcycler_df <- load_phenocycler_summary_files()
#' }
#' @export
load_phenocycler_summary_files <- function(data_dir = NULL) {
    pcycler_dt <- .find_and_load_phenocycler_summary_files(data_dir = data_dir)

    meta_df <- load_meta_data()


    pcycler_dt = pcycler_dt %>% mutate(sample_id = `Image Tag`) %>%
        mutate(sample_id = sub("\\..+", "", sample_id)) %>%
        mutate(sample_id = gsub("-", "_", sample_id)) %>%
        mutate(sample_id = sub("DEB", "D_EB_", sample_id)) %>%
        mutate(sample_id = sub("CTEBV", "CTEBV_", sample_id)) %>%
        mutate(probe_control = ifelse(grepl("Control", sample_id), "negative_probe", "")) %>%
        mutate(sample_id = sub("CellPelletSlide_Control_Scan1_Phenocycler_", "", sample_id)) %>%
        mutate(sample_id = sub("CellPelletSlide_Test_Scan1_Phenocycler_", "", sample_id)) %>%
        mutate(sample_id = sub("_Scan.+", "", sample_id))


    setdiff(pcycler_dt$sample_id, meta_df$sample_id)


    pcycler_dt <- merge(pcycler_dt, meta_df, all.x = TRUE, by = "sample_id")
    pcycler_dt <- dplyr::mutate(
        pcycler_dt,
        EBER_status = ifelse(is.na(.data$EBER_status), "need info", .data$EBER_status)
    )
    pcycler_dt$assay = EBV_ASSAY_TYPES$Phenocycler
    pcycler_dt$project_name = assay_to_project_name[EBV_ASSAY_TYPES$Phenocycler]
    pcycler_dt
}
