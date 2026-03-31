.validate_CellDataStore <- function(object) {
  if (!all(c("sample_id", "n_cells_original") %in% colnames(object@original_cells_per_sample_id))) {
    return("`original_cells_per_sample_id` must contain columns `sample_id` and `n_cells_original`.")
  }
  TRUE
}

#' Cell Data Store Container
#'
#' S4 class that stores the full loaded cell table alongside a mutable selected
#' view used for iterative filtering and transformation.
#'
#' @slot full_df Full cell-level data frame.
#' @slot selected_df Selected cell-level data frame used by dplyr operations.
#' @slot original_cells_per_sample_id Data frame with original per-sample cell
#'   counts at object creation time.
#'
#' @exportClass CellDataStore
methods::setClass(
  "CellDataStore",
  slots = c(
    full_df = "data.frame",
    selected_df = "data.frame",
    original_cells_per_sample_id = "data.frame"
  ),
  validity = .validate_CellDataStore
)

.count_original_cells <- function(cell_df) {
  if (!nrow(cell_df)) {
    return(data.frame(sample_id = character(0), n_cells_original = integer(0)))
  }

  stopifnot("sample_id" %in% colnames(cell_df))

  out <- as.data.frame(table(cell_df$sample_id), stringsAsFactors = FALSE)
  colnames(out) <- c("sample_id", "n_cells_original")
  out$n_cells_original <- as.integer(out$n_cells_original)
  out
}

#' Construct a CellDataStore Object
#'
#' Creates a [CellDataStore-class] where `selected_df` is initialized to
#' `full_df`, and original per-sample cell counts are stored in
#' `original_cells_per_sample_id`.
#'
#' @param cell_df Cell-level data frame with a `sample_id` column.
#'
#' @return A validated [CellDataStore-class] object.
#' @examples
#' \dontrun{
#' q <- CellQuery()
#' cell_df <- load_query_cell_data(q)
#' cds <- new_cell_data_store(cell_df)
#' }
#' @export
new_cell_data_store <- function(cell_df) {
  stopifnot(is.data.frame(cell_df))
  stopifnot("sample_id" %in% colnames(cell_df))

  obj <- methods::new(
    "CellDataStore",
    full_df = cell_df,
    selected_df = cell_df,
    original_cells_per_sample_id = .count_original_cells(cell_df)
  )
  methods::validObject(obj)
  obj
}

#' Load Query Cell Data As CellDataStore
#'
#' Convenience wrapper around [load_query_cell_data()] + [new_cell_data_store()].
#'
#' @param object A [CellQueryInfo-class] object.
#'
#' @return A [CellDataStore-class] object.
#' @examples
#' \dontrun{
#' q <- CellQuery()
#' cds <- load_query_cell_data_store(q)
#' }
#' @export
load_query_cell_data_store <- function(object) {
  cell_df <- load_query_cell_data(object)
  new_cell_data_store(cell_df)
}

#' Get Full Cell Data
#'
#' @param object A [CellDataStore-class] object.
#'
#' @return The full cell-level data frame.
#' @examples
#' \dontrun{
#' full_df <- get_full_cell_data(cds)
#' }
#' @export
get_full_cell_data <- function(object) {
  stopifnot(methods::is(object, "CellDataStore"))
  object@full_df
}

#' Get Selected Cell Data
#'
#' @param object A [CellDataStore-class] object.
#'
#' @return The selected cell-level data frame.
#' @examples
#' \dontrun{
#' sel_df <- get_selected_cell_data(cds)
#' }
#' @export
get_selected_cell_data <- function(object) {
  stopifnot(methods::is(object, "CellDataStore"))
  object@selected_df
}

#' Get Original Cell Counts Per Sample
#'
#' @param object A [CellDataStore-class] object.
#'
#' @return Data frame with `sample_id` and `n_cells_original` columns.
#' @examples
#' \dontrun{
#' get_original_cells_per_sample_id(cds)
#' }
#' @export
get_original_cells_per_sample_id <- function(object) {
  stopifnot(methods::is(object, "CellDataStore"))
  object@original_cells_per_sample_id
}

#' Restore Selected Data To Full Data
#'
#' Resets `selected_df` to the complete `full_df`.
#'
#' @param object A [CellDataStore-class] object.
#'
#' @return Updated [CellDataStore-class] object.
#' @examples
#' \dontrun{
#' cds <- restore_selected_cell_data(cds)
#' }
#' @export
restore_selected_cell_data <- function(object) {
  stopifnot(methods::is(object, "CellDataStore"))
  object@selected_df <- object@full_df
  methods::validObject(object)
  object
}

#' Replace Full Data With Selected Data
#'
#' Commits the currently selected view by replacing `full_df` with
#' `selected_df`.
#'
#' @param object A [CellDataStore-class] object.
#'
#' @return Updated [CellDataStore-class] object.
#' @examples
#' \dontrun{
#' cds <- replace_full_with_selected_cell_data(cds)
#' }
#' @export
replace_full_with_selected_cell_data <- function(object) {
  stopifnot(methods::is(object, "CellDataStore"))
  object@full_df <- object@selected_df
  methods::validObject(object)
  object
}

#' Display a CellDataStore Summary
#'
#' @param object A [CellDataStore-class] object.
#'
#' @return Invisibly returns `object`.
#' @export
methods::setMethod("show", "CellDataStore", function(object) {
  cat("CellDataStore\n")
  cat("  full rows:", nrow(object@full_df), "\n")
  cat("  selected rows:", nrow(object@selected_df), "\n")
  cat("  original sample count rows:", nrow(object@original_cells_per_sample_id), "\n")
  invisible(object)
})

#' @importFrom dplyr filter
#' @export
filter.CellDataStore <- function(.data, ...) {
  .data@selected_df <- dplyr::filter(.data@selected_df, ...)
  .data
}

#' @importFrom dplyr mutate
#' @export
mutate.CellDataStore <- function(.data, ...) {
  .data@selected_df <- dplyr::mutate(.data@selected_df, ...)
  .data
}

#' @importFrom dplyr select
#' @export
select.CellDataStore <- function(.data, ...) {
  .data@selected_df <- dplyr::select(.data@selected_df, ...)
  .data
}

#' @importFrom dplyr arrange
#' @export
arrange.CellDataStore <- function(.data, ...) {
  .data@selected_df <- dplyr::arrange(.data@selected_df, ...)
  .data
}

#' @importFrom dplyr slice
#' @export
slice.CellDataStore <- function(.data, ..., .by = NULL, .preserve = FALSE) {
  .data@selected_df <- dplyr::slice(.data@selected_df, ..., .by = .by, .preserve = .preserve)
  .data
}

#' @importFrom dplyr distinct
#' @export
distinct.CellDataStore <- function(.data, ..., .keep_all = FALSE) {
  .data@selected_df <- dplyr::distinct(.data@selected_df, ..., .keep_all = .keep_all)
  .data
}
