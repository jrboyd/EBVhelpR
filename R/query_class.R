
.validate_CellQueryInfo = function(object) {
    if (length(object@assay_type) != 1L || is.na(object@assay_type) || object@assay_type == "") {
        return("`assay_type` must be a non-empty character scalar.")
    }
    if (!is.character(object@selected_sample_ids)) {
        return("`selected_sample_ids` must be a character vector.")
    }
    all_sample_ids <- .get_sample_ids(object@all_cell_files_df)
    if (length(all_sample_ids)) {
        is_valid <- object@selected_sample_ids %in% all_sample_ids
        if (!all(is_valid)) {
            return("All `selected_sample_ids` must be present in `all_cell_files_df`.")
        }
    }
    TRUE
}

#' Cell Query Information Container
#'
#' S4 class for storing data used during cell-query image retrieval workflows.
#'
#' @slot all_cell_files_df Data frame of all available cell files.
#' @slot meta_data_df Data frame returned by [load_meta_data()].
#' @slot selected_sample_ids Character vector of selected sample ids.
#' @slot tiff_paths_df Data frame mapping samples to TIFF image paths.
#' @slot assay_type Character scalar indicating assay type (for example,
#'   `"RNAScope_4plex"`).
#'
#' @exportClass CellQueryInfo
methods::setClass(
    "CellQueryInfo",
    slots = c(
        summary_df = "data.frame",
        all_cell_files_df = "data.frame",
        meta_data_df = "data.frame",
        selected_sample_ids = "character",
        selected_unique_ids = "character",
        tiff_paths_df = "data.frame",
        assay_type = "character"
    ),
    validity = .validate_CellQueryInfo
)

.df_prep = function(df){
    stopifnot("sample_id" %in% colnames(df))
    stopifnot("probe_control" %in% colnames(df))
    df = df %>% dplyr::mutate(unique_id =
                                  ifelse(probe_control == "",
                                         sample_id,
                                         paste(sample_id, probe_control)))
    df
}

#' Construct a CellQueryInfo Object
#'
#' Builds assay-specific summary, cell file, and TIFF path tables, and
#' initializes selection to all sample ids in the cell-file table.
#'
#' @param assay_type Character scalar assay type. If `NULL`, defaults to
#'   `EBV_ASSAY_TYPES$RNAScope_4plex`.
#'
#' @return A validated [CellQueryInfo-class] object.
#' @examples
#' q <- CellQuery()
#' q2 <- CellQuery(EBV_ASSAY_TYPES$Phenocycler)
#' @export
CellQuery <- function(
        assay_type = NULL
) {
    if(is.null(assay_type)){
        assay_type = EBV_ASSAY_TYPES$RNAScope_4plex
        message("Defaulting to assay type ", assay_type, ".")
        message("Select valid assay types with EBV_ASSAY_TYPES, i.e. EBV_ASSAY_TYPES$RNAScope_4plex")
    }
    stopifnot(assay_type %in% EBV_ASSAY_TYPES)
    if(assay_type == EBV_ASSAY_TYPES$Phenocycler){
        summary_df = load_phenocycler_summary_files()
    }else if(assay_type %in% c(EBV_ASSAY_TYPES$RNAScope_4plex, EBV_ASSAY_TYPES$`RNAScope_3plex+IF`)){
        summary_df = load_rnascope_summary_files()
        summary_df = dplyr::filter(summary_df, assay == assay_type)
    }else{
        stop("Unrecognized assay type, see EBV_ASSAY_TYPES")
    }
    all_cell_files_df = load_cell_source_files()
    all_cell_files_df = dplyr::filter(all_cell_files_df, assay == assay_type)

    tiff_paths_df = get_tiff_file_path_df()
    tiff_paths_df = dplyr::filter(tiff_paths_df, assay == assay_type)

    meta_data_df <- load_meta_data()

    .warn_undocumented_sample_ids(
        meta_data_df = meta_data_df,
        dfs = list(
            summary_df = summary_df,
            all_cell_files_df = all_cell_files_df,
            tiff_paths_df = tiff_paths_df
        )
    )

    selected_sample_ids = all_cell_files_df$sample_id %>% unique

    summary_df = .df_prep(summary_df)
    all_cell_files_df = .df_prep(all_cell_files_df)
    tiff_paths_df = .df_prep(tiff_paths_df)

    selected_unique_ids = all_cell_files_df$sample_id %>% unique

    obj <- methods::new(
        "CellQueryInfo",
        summary_df = summary_df,
        all_cell_files_df = all_cell_files_df,
        meta_data_df = meta_data_df,
        selected_sample_ids = selected_sample_ids,
        selected_unique_ids = selected_unique_ids,
        tiff_paths_df = tiff_paths_df,
        assay_type = assay_type
    )
    methods::validObject(obj)
    obj
}

.resolve_col_name <- function(df, candidates) {
    idx <- candidates[candidates %in% colnames(df)]
    if (!length(idx)) {
        return(NULL)
    }
    idx[[1]]
}

.get_sample_control_ids = function(df){
    sample_col = "sample_id"
    control_col = "probe_control"
    if (is.null(sample_col) || !nrow(df)) {
        return(character(0))
    }
    unique(
        paste(
            as.character(stats::na.omit(df[[sample_col]])),
            as.character(stats::na.omit(df[[control_col]]))
        )
    )
}

.get_sample_ids <- function(df, sample_col = "sample_id") {

    if (is.null(sample_col) || !nrow(df)) {
        return(character(0))
    }
    unique(as.character(stats::na.omit(df[[sample_col]])))
}

.filter_by_selected_sample_ids <- function(df, selected_sample_ids, sample_col = "sample_id") {
    if (!nrow(df)) {
        return(df)
    }
    if (!length(selected_sample_ids)) {
        return(df[0, , drop = FALSE])
    }
    df[df[[sample_col]] %in% selected_sample_ids, , drop = FALSE]
}

.warn_undocumented_sample_ids <- function(meta_data_df, dfs) {
    documented_ids <- .get_sample_ids(meta_data_df)
    if (!length(documented_ids)) {
        warning("No metadata sample_id values were found; unable to check for undocumented sample ids.")
        return(invisible(NULL))
    }

    for (df_name in names(dfs)) {
        df_sample_ids <- .get_sample_ids(dfs[[df_name]])
        undocumented <- setdiff(df_sample_ids, documented_ids)
        if (length(undocumented)) {
            warning(
                "Undocumented sample_id values found in ", df_name, ": ",
                paste(undocumented, collapse = ", ")
            )
        }
    }

    invisible(NULL)
}

.limit_query_to_sample_ids <- function(object, sample_ids, unique_ids) {
    sample_ids <- unique(as.character(sample_ids))

    object@meta_data_df <- .filter_by_selected_sample_ids(object@meta_data_df, sample_ids)

    object@summary_df <- .filter_by_selected_sample_ids(
        object@summary_df,
        unique_ids, "unique_id"
    )
    object@all_cell_files_df <- .filter_by_selected_sample_ids(
        object@all_cell_files_df,
        unique_ids, "unique_id"
    )
    object@tiff_paths_df <- .filter_by_selected_sample_ids(
        object@tiff_paths_df,
        unique_ids, "unique_id"
    )

    object@selected_sample_ids <- intersect(object@selected_sample_ids, sample_ids)
    object@selected_sample_ids <- unique(as.character(object@selected_sample_ids))

    object@selected_unique_ids <- intersect(object@selected_unique_ids, unique_ids)
    object@selected_unique_ids <- unique(as.character(object@selected_unique_ids))

    methods::validObject(object)
    object
}

#' Filter Query Data To Cell-File Sample IDs
#'
#' Restricts all stored data frames in a [CellQueryInfo-class] object to
#' sample ids present in `all_cell_files_df`.
#'
#' @param object A [CellQueryInfo-class] object.
#'
#' @return Updated [CellQueryInfo-class] object.
#' @export
filter_query_to_all_cell_file_samples <- function(object) {
    stopifnot(methods::is(object, "CellQueryInfo"))
    .limit_query_to_sample_ids(
        object,
        .get_sample_ids(object@all_cell_files_df),
        .get_sample_ids(object@all_cell_files_df, sample_col = "unique_id")
    )
}

#' Filter Query Data To TIFF-Path Sample IDs
#'
#' Restricts all stored data frames in a [CellQueryInfo-class] object to
#' sample ids present in `tiff_paths_df`.
#'
#' @param object A [CellQueryInfo-class] object.
#'
#' @return Updated [CellQueryInfo-class] object.
#' @export
filter_query_to_tiff_path_samples <- function(object) {
    stopifnot(methods::is(object, "CellQueryInfo"))
    .limit_query_to_sample_ids(
        object,
        .get_sample_ids(object@tiff_paths_df),
        .get_sample_ids(object@tiff_paths_df, sample_col = "unique_id")
    )
}

#' Get Summary Rows From CellQueryInfo
#'
#' @param object A [CellQueryInfo-class] object.
#' @param selected_only Logical; if `TRUE`, return only selected samples.
#'
#' @return Data frame of summary rows.
#' @examples
#' q <- CellQuery()
#' get_query_summary_df(q)
#' @export
get_query_summary_df <- function(object, selected_only = TRUE) {
    stopifnot(methods::is(object, "CellQueryInfo"))
    if (!selected_only) {
        return(object@summary_df)
    }
    .filter_by_selected_sample_ids(object@summary_df, object@selected_sample_ids)
}

#' Get Cell File Rows From CellQueryInfo
#'
#' @param object A [CellQueryInfo-class] object.
#' @param selected_only Logical; if `TRUE`, return only selected samples.
#'
#' @return Data frame of cell file rows.
#' @examples
#' q <- CellQuery()
#' get_query_cell_files_df(q)
#' @export
get_query_cell_files_df <- function(object, selected_only = TRUE) {
    stopifnot(methods::is(object, "CellQueryInfo"))
    if (!selected_only) {
        return(object@all_cell_files_df)
    }
    .filter_by_selected_sample_ids(object@all_cell_files_df, object@selected_unique_ids, sample_col = "unique_id")
}

#' Get TIFF Path Rows From CellQueryInfo
#'
#' @param object A [CellQueryInfo-class] object.
#' @param selected_only Logical; if `TRUE`, return only selected samples.
#'
#' @return Data frame of TIFF path rows.
#' @examples
#' q <- CellQuery()
#' get_query_tiff_paths_df(q)
#' @export
get_query_tiff_paths_df <- function(object, selected_only = TRUE) {
    stopifnot(methods::is(object, "CellQueryInfo"))
    if (!selected_only) {
        return(object@tiff_paths_df)
    }
    .filter_by_selected_sample_ids(object@tiff_paths_df, object@selected_sample_ids)
}

#' Set Selected Sample IDs
#'
#' @param object A [CellQueryInfo-class] object.
#' @param selected_sample_ids Character vector of selected sample ids.
#'
#' @return Updated [CellQueryInfo-class] object.
#' @examples
#' q <- CellQuery()
#' q <- set_selected_sample_ids(q, head(get_query_summary_df(q)$sample_id, 2))
#' @export
set_selected_sample_ids <- function(object, selected_sample_ids) {
    stopifnot(methods::is(object, "CellQueryInfo"))
    object@selected_sample_ids <- unique(as.character(selected_sample_ids))
    methods::validObject(object)
    object
}

#' Set Selected Sample IDs
#'
#' @param object A [CellQueryInfo-class] object.
#' @param selected_sample_ids Character vector of selected sample ids.
#'
#' @return Updated [CellQueryInfo-class] object.
#' @examples
#' q <- CellQuery()
#' q <- set_selected_sample_ids(q, head(get_query_summary_df(q)$sample_id, 2))
#' @export
set_selected_unique_ids <- function(object, selected_unique_ids) {
    stopifnot(methods::is(object, "CellQueryInfo"))
    object@selected_unique_ids <- unique(as.character(selected_unique_ids))
    methods::validObject(object)
    object
}

.count_samples <- function(df) {
    sample_col = "sample_id"
    if (is.null(sample_col)) {
        return(0L)
    }
    length(unique(df[[sample_col]]))
}

.count_cells <- function(df) {
    if (!nrow(df)) {
        return(0L)
    }
    cell_count_col <- .resolve_col_name(df, c("n_cells", "cell_count", "N_cells", "n"))
    if (!is.null(cell_count_col) && is.numeric(df[[cell_count_col]])) {
        return(as.integer(sum(df[[cell_count_col]], na.rm = TRUE)))
    }
    nrow(df)
}

.count_tiffs_for_samples <- function(tiff_df, selected_sample_ids) {
    if (!length(selected_sample_ids) || !nrow(tiff_df)) {
        return(0L)
    }
    sample_col <- .resolve_col_name(tiff_df, c("sample_id", "SampleID", "Sample", "name"))
    tiff_col <- .resolve_col_name(tiff_df, c("tiff_file", "tiff_path", "file", "path"))
    if (is.null(sample_col)) {
        return(0L)
    }

    sel <- tiff_df[[sample_col]] %in% selected_sample_ids
    if (!any(sel)) {
        return(0L)
    }

    if (is.null(tiff_col)) {
        return(length(unique(tiff_df[[sample_col]][sel])))
    }

    has_path <- !is.na(tiff_df[[tiff_col]]) & nzchar(tiff_df[[tiff_col]])
    length(unique(tiff_df[[tiff_col]][sel & has_path]))
}

.report_sample_overlaps = function(cq){
    samp_ids.l = list(
        documented = cq@meta_data_df$sample_id,
        Summarized = cq@summary_df$sample_id,
        `Cell Data` = cq@all_cell_files_df$sample_id,
        `Tiffs Available` = cq@tiff_paths_df$sample_id
    )
    uniq_ids.l = list(
        Summarized = cq@summary_df$unique_id,
        `Cell Data` = cq@all_cell_files_df$unique_id,
        `Tiffs Available` = cq@tiff_paths_df$unique_id
    )

    samp_ids.l = lapply(samp_ids.l, unique)
    uniq_ids.l = lapply(uniq_ids.l, unique)

    assayed_ids = unlist(samp_ids.l[-1]) %>% unique
    perc_assayed = round(100*length(assayed_ids) / length(samp_ids.l$documented), 1)

    .report_overlap = function(x, other){
        paste0(length(x), " (", round(100*length(x)/length(other), 1), "%)")
    }

    unique_ids = unlist(uniq_ids.l) %>% unique
    slot_reports = sapply(uniq_ids.l, function(x){
        .report_overlap(x, unique_ids)
    })

    message(length(samp_ids.l$documented), " total documented samples")
    message(.report_overlap(assayed_ids, samp_ids.l$documented), " have assay data")
    message(length(unique_ids), " total samples, counting controls")
    message("of assayed samples:")
    for(name in names(slot_reports)){
        pad_len = max(nchar(slot_reports)) + 2
        pad = paste(rep(" ", pad_len - nchar(slot_reports[name])), collapse = "")
        message("  ", slot_reports[name],  pad, ":", name)
    }
}

#' Display a CellQueryInfo Summary
#'
#' Prints assay name and key counts for total versus selected samples, cells,
#' and available TIFF files.
#'
#' @param object A [CellQueryInfo-class] object.
#'
#' @return Invisibly returns `object`.
#' @examples
#' q <- CellQuery()
#' show(q)
#' @export
methods::setMethod("show", "CellQueryInfo", function(object) {
    # all_df <- get_query_cell_files_df(object, selected_only = FALSE)
    # sel_df <- get_query_cell_files_df(object, selected_only = TRUE)
    # tiff_df <- get_query_tiff_paths_df(object, selected_only = FALSE)
    # assay <- object@assay_type
    #
    # total_samples <- .count_samples(all_df)
    # selected_samples <- length(unique(object@selected_sample_ids))
    # total_cells <- .count_cells(all_df)
    # selected_cells <- .count_cells(sel_df)
    # selected_tiffs <- .count_tiffs_for_samples(tiff_df, object@selected_sample_ids)
    #
    # cat("CellQueryInfo\n")
    # cat("  assay:", assay, "\n")
    # cat("  total samples available:", total_samples, "\n")
    # cat("  selected samples:", selected_samples, "\n")
    # cat("  total cells:", total_cells, "\n")
    # cat("  selected cells:", selected_cells, "\n")
    # cat("  tiffs for selected samples:", selected_tiffs, "\n")
    # invisible(object)
    .report_sample_overlaps(object)
})


