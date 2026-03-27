#' Cell Query Information Container
#'
#' S4 class for storing data used during cell-query image retrieval workflows.
#'
#' @slot all_cell_files_df Data frame of all available cell files.
#' @slot selected_cell_files_df Data frame of cell files selected for a query.
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
        selected_cell_files_df = "data.frame",
        tiff_paths_df = "data.frame",
        assay_type = "character"
    ),
    validity = function(object) {
        if (length(object@assay_type) != 1L || is.na(object@assay_type) || object@assay_type == "") {
            return("`assay_type` must be a non-empty character scalar.")
        }
        TRUE
    }
)

#' Construct a CellQueryInfo Object
#'
#' @param all_cell_files_df Data frame of all available cell files.
#' @param selected_cell_files_df Data frame of selected cell files.
#' @param tiff_paths_df Data frame with TIFF path information.
#' @param assay_type Character scalar indicating assay type.
#'
#' @return A validated [CellQueryInfo-class] object.
#' @export
#'
#' @examples
#' assay_type = EBV_ASSAY_TYPES$phenocycler
#' all_cell_files_df = NULL
#' selected_cell_files_df = all_cell_files_df
#' tiff_paths_df = NULL
#' get_tiff_file_path_df = EBVhelpR:::get_tiff_file_path_df
new_cell_query_info <- function(
        assay_type = NULL,
        summary_df = NULL,
        all_cell_files_df = NULL,
        selected_cell_files_df = all_cell_files_df,
        tiff_paths_df = NULL
) {
    if(is.null(assay_type)){
       assay_type = EBV_ASSAY_TYPES$rnascope_4plex
       message("Defaulting to assay type ", assay_type, ".")
       message("Select valid assay types with EBV_ASSAY_TYPES, i.e. EBV_ASSAY_TYPES$rnascope_4plex.")
    }
    stopifnot(assay_type %in% EBV_ASSAY_TYPES)
    if(assay_type == EBV_ASSAY_TYPES$phenocycler){
        summary_df = load_phenocycler_summary_files()
        all_cell_files_df = load_cell_source_files()
        all_cell_files_df = dplyr::filter(all_cell_files_df, assay == assay_type)
        selected_cell_files_df = all_cell_files_df
        tiff_info = get_tiff_file_path_df()
        tiff_info = dplyr::filter(tiff_info, assay == assay_type)
    }else if(assay_type %in% c(EBV_ASSAY_TYPES$rnascope_4plex, EBV_ASSAY_TYPES$`rnascope_3plex+IF`)){
        summary_df = load_rnascope_summary_files()
        summary_df = dplyr::filter(summary_df, assay == assay_type)
        all_cell_files_df = load_cell_source_files()
        all_cell_files_df = dplyr::filter(all_cell_files_df, assay == assay_type)
        selected_cell_files_df = all_cell_files_df
        tiff_info = get_tiff_file_path_df()
        tiff_info = dplyr::filter(tiff_info, assay == assay_type)
    }else{
        stop("Unrecognized assay type, see EBV_ASSAY_TYPES")
    }
    obj <- methods::new(
        "CellQueryInfo",
        summary_df = summary_df,
        all_cell_files_df = all_cell_files_df,
        selected_cell_files_df = selected_cell_files_df,
        tiff_paths_df = tiff_paths_df,
        assay_type = assay_type
    )
    methods::validObject(obj)
    obj
}


