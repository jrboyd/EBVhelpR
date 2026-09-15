#' Phenocycler gate column renaming
#'
#' The Phenocycler export names gate columns with `+` and `-` — `CD4+TCells` and
#' `CD4-Tcells` are a population and its complement. `make.names()` maps **both**
#' characters to `.`, so the pair collapses to `CD4.TCells`/`CD4.Tcells`,
#' distinguished only by the case of the `T`, and `CD8+Tcells`/`CD8-Tcells`
#' become `CD8.Tcells` and `CD8.Tcells.1`. A single case typo then analyses the
#' opposite cells.
#'
#' Renaming happens at import, which is the last point where `+` and `-` are
#' still present: `readr::read_csv()` preserves them, but any downstream
#' `as.data.frame()`, `read.csv()` or other `check.names = TRUE` path does not.
#'
#' `PDL1+Bcells` is gated on PDL1 and **Pax5**, not CD20; the new name records
#' that.
#'
#' @noRd
.phenocycler_gate_renames <- function() {
    base_map <- c(
        "CD4+TCells"  = "CD4_pos_Tcells",
        "CD4-Tcells"  = "CD4_neg_Tcells",
        "CD8+Tcells"  = "CD8_pos_Tcells",
        "CD8-Tcells"  = "CD8_neg_Tcells",
        "PD1+Tcells"  = "PD1_pos_Tcells",
        "PDL1+Bcells" = "PDL1_pos_Pax5_Bcells"
    )
    # Each gate appears as a count column and as a percentage column.
    counts <- stats::setNames(
        paste0(base_map, " Cells"),
        paste0(names(base_map), " Cells")
    )
    percents <- stats::setNames(
        paste0("% ", base_map, " Positive Cells"),
        paste0("% ", names(base_map), " Positive Cells")
    )
    c(counts, percents)
}

#' Rename `+`/`-` Phenocycler gate columns to unambiguous names
#'
#' @param dt Phenocycler summary data frame.
#' @return `dt` with gate columns renamed.
#' @noRd
.rename_phenocycler_gate_columns <- function(dt) {
    renames <- .phenocycler_gate_renames()
    cn <- colnames(dt)

    missing <- setdiff(names(renames), cn)
    if (length(missing)) {
        stop(
            "Phenocycler summaries are missing expected gate column(s): ",
            paste(missing, collapse = ", "),
            ". The export format may have changed.",
            call. = FALSE
        )
    }

    idx <- match(names(renames), cn)
    cn[idx] <- unname(renames)
    colnames(dt) <- cn
    dt
}

#' Assert the CD4 and CD8 gates are complementary partitions of CD3+ cells
#'
#' `CD4_pos + CD4_neg` and `CD8_pos + CD8_neg` are each the CD3+ count, so they
#' must be equal. The summary export carries no CD3 column of its own, which is
#' why the check is written as the equality of the two sums. This holds exactly
#' across the current data; a failure means the gating definitions changed.
#'
#' @param dt Phenocycler summary data frame with renamed gate columns.
#' @return `invisible(NULL)`, called for the error.
#' @noRd
.assert_phenocycler_gate_complement <- function(dt) {
    cd4 <- dt[["CD4_pos_Tcells Cells"]] + dt[["CD4_neg_Tcells Cells"]]
    cd8 <- dt[["CD8_pos_Tcells Cells"]] + dt[["CD8_neg_Tcells Cells"]]

    bad <- which(!is.na(cd4) & !is.na(cd8) & cd4 != cd8)
    if (length(bad)) {
        tags <- as.character(dt[["Image Tag"]])[bad]
        stop(
            "CD4 and CD8 gates do not partition the same CD3+ population in ",
            length(bad), " row(s): ",
            paste(utils::head(tags, 10), collapse = ", "),
            if (length(bad) > 10) ", ..." else "",
            ". Expected `CD4_pos + CD4_neg == CD8_pos + CD8_neg`.",
            call. = FALSE
        )
    }
    invisible(NULL)
}

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

    dt <- .rename_phenocycler_gate_columns(dt)
    .assert_phenocycler_gate_complement(dt)

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
