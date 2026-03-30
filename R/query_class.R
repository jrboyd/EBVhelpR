
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
        selected_sample_ids = "character",
        tiff_paths_df = "data.frame",
        assay_type = "character"
    ),
    validity = .validate_CellQueryInfo
)

#' Construct a CellQueryInfo Object
#'
#' Builds assay-specific summary, cell file, and TIFF path tables; keeps the
#' intersection of available sample ids across those sources; and initializes
#' selection to all remaining sample ids.
#'
#' @param assay_type Character scalar assay type. If `NULL`, defaults to
#'   `EBV_ASSAY_TYPES$rnascope_4plex`.
#'
#' @return A validated [CellQueryInfo-class] object.
#' @examples
#' q <- CellQuery()
#' q2 <- CellQuery(EBV_ASSAY_TYPES$phenocycler)
#' @export
CellQuery <- function(
        assay_type = NULL
) {
    if(is.null(assay_type)){
       assay_type = EBV_ASSAY_TYPES$rnascope_4plex
       message("Defaulting to assay type ", assay_type, ".")
       message("Select valid assay types with EBV_ASSAY_TYPES, i.e. EBV_ASSAY_TYPES$rnascope_4plex.")
    }
    stopifnot(assay_type %in% EBV_ASSAY_TYPES)
    if(assay_type == EBV_ASSAY_TYPES$phenocycler){
        summary_df = load_phenocycler_summary_files()
    }else if(assay_type %in% c(EBV_ASSAY_TYPES$rnascope_4plex, EBV_ASSAY_TYPES$`rnascope_3plex+IF`)){
        summary_df = load_rnascope_summary_files()
        summary_df = dplyr::filter(summary_df, assay == assay_type)
    }else{
        stop("Unrecognized assay type, see EBV_ASSAY_TYPES")
    }
    all_cell_files_df = load_cell_source_files()
    all_cell_files_df = dplyr::filter(all_cell_files_df, assay == assay_type)

    tiff_paths_df = get_tiff_file_path_df()
    tiff_paths_df = dplyr::filter(tiff_paths_df, assay == assay_type)

    common_ids = intersect(summary_df$sample_id, all_cell_files_df$sample_id)
    common_ids = intersect(common_ids, tiff_paths_df$sample_id)

    summary_drop = setdiff(summary_df$sample_id, common_ids)
    cell_drop = setdiff(all_cell_files_df$sample_id, common_ids)
    tiff_drop = setdiff(tiff_paths_df$sample_id, common_ids)

    if(length(summary_drop) > 0){
        warning("dropping ids from summary:\n", paste(summary_drop, collapse = ", "))
    }
    if(length(cell_drop) > 0){
        warning("dropping ids from cell info:\n", paste(cell_drop, collapse = ", "))
    }
    if(length(tiff_drop) > 0){
        warning("dropping ids from tiff paths:\n", paste(tiff_drop, collapse = ", "))
    }

    summary_df = summary_df %>%
        dplyr::filter(sample_id %in% common_ids)
    all_cell_files_df = all_cell_files_df %>%
        dplyr::filter(sample_id %in% common_ids)
    tiff_paths_df = tiff_paths_df %>%
        dplyr::filter(sample_id %in% common_ids)

    selected_sample_ids = all_cell_files_df$sample_id %>% unique

    obj <- methods::new(
        "CellQueryInfo",
        summary_df = summary_df,
        all_cell_files_df = all_cell_files_df,
        selected_sample_ids = selected_sample_ids,
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

.get_sample_ids <- function(df) {
    sample_col = "sample_id"
    if (is.null(sample_col) || !nrow(df)) {
        return(character(0))
    }
    unique(as.character(stats::na.omit(df[[sample_col]])))
}

.filter_by_selected_sample_ids <- function(df, selected_sample_ids) {
    if (!nrow(df)) {
        return(df)
    }
    sample_col = "sample_id"
    if (!length(selected_sample_ids)) {
        return(df[0, , drop = FALSE])
    }
    df[df[[sample_col]] %in% selected_sample_ids, , drop = FALSE]
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
    .filter_by_selected_sample_ids(object@all_cell_files_df, object@selected_sample_ids)
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
    all_df <- get_query_cell_files_df(object, selected_only = FALSE)
    sel_df <- get_query_cell_files_df(object, selected_only = TRUE)
    tiff_df <- get_query_tiff_paths_df(object, selected_only = FALSE)
    assay <- object@assay_type

    total_samples <- .count_samples(all_df)
    selected_samples <- length(unique(object@selected_sample_ids))
    total_cells <- .count_cells(all_df)
    selected_cells <- .count_cells(sel_df)
    selected_tiffs <- .count_tiffs_for_samples(tiff_df, object@selected_sample_ids)

    cat("CellQueryInfo\n")
    cat("  assay:", assay, "\n")
    cat("  total samples available:", total_samples, "\n")
    cat("  selected samples:", selected_samples, "\n")
    cat("  total cells:", total_cells, "\n")
    cat("  selected cells:", selected_cells, "\n")
    cat("  tiffs for selected samples:", selected_tiffs, "\n")
    invisible(object)
})


