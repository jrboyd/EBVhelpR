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
new_cell_query_info <- function(
        assay_type,
        all_cell_files_df = NULL,
        selected_cell_files_df = all_cell_files_df,
        tiff_paths_df = NULL,

) {
    stopifnot(assay_type %in% EBV_ASSAY_TYPES)
    if(assay_type == EBV_ASSAY_TYPES$phenocycler){
        all_cell_files_df = load_phenocycler_summary_files()
        cell_info = load_cell_source_files()
        cell_info = subset(cell_info, assay == assay_type)
    }else if(assay_type == EBV_ASSAY_TYPES$rnascope_4plex){

    }else if(assay_type == EBV_ASSAY_TYPES$`rnascope_3plex+IF`){

    }
    obj <- methods::new(
        "CellQueryInfo",
        all_cell_files_df = all_cell_files_df,
        selected_cell_files_df = selected_cell_files_df,
        tiff_paths_df = tiff_paths_df,
        assay_type = assay_type
    )
    methods::validObject(obj)
    obj
}

EBV_ASSAY_TYPES = list(
    phenocycler = "phenocycler",
    rnascope_4plex = "rnascope_4plex",
    "rnascope_3plex+IF" = "rnascope_3plex+IF"
)

EBV_DIR_NAMES = list(
    phenocycler = "Phenocycler",
    rnascope_4plex = "RNAScopeRound1",
    "rnascope_3plex+IF" = "RNAScopeRound2"
)

.get_valid_project_names = function(){
    unlist(EBV_DIR_NAMES[unlist(EBV_ASSAY_TYPES)])
}

get_tiff_file_path_df = function(){
    win_dir = "Z:/FUSION DATA/AshleyVolaric"
    lin_dir = "/netfiles/volaric_research/DLBCL_EBV_detection/image_files"
    tiff_dir = win_dir
    if(!dir.exists(tiff_dir)){
        tiff_dir = lin_dir
    }
    if(!dir.exists(tiff_dir)){
        stop("Could not locate TIFF root directory.")
    }
    dir_names = dir(tiff_dir)
    stopifnot(all(.get_valid_project_names() %in% dir_names))
    dir_names = .get_valid_project_names()
    names(dir_names) = dir_names
    tiff_files.by_project = lapply(dir_names, function(d){
        files = dir(file.path(tiff_dir, d), recursive = TRUE, pattern = "tiff?$", full.names = TRUE)
        data.frame(tiff_file = files)
    })
    tiff_df = bind_rows(tiff_files.by_project, .id = "assay")

    tiff_df = tiff_df %>% mutate(name = basename(tiff_file))
    tiff_df = tiff_df %>%
        mutate(name = sub("\\..+", "", name)) %>%
        mutate(name = sub("^[0-9]{3}_", "", name)) %>%
        mutate(name = sub("_ ?croo?p?ped$", "", name)) %>%
        mutate(name = gsub("-", "_", name))

    tiff_df = tiff_df %>%
        mutate(name = ifelse(
            grepl("D_EB.+D_EB", name),
            sub("_D", " D", name),
            name
        ))
    tiff_df = tiff_df %>%
        mutate(name = ifelse(
            grepl("CTEBV.+CTEBV", name),
            sub("_CTEBV", " CTEBV", name),
            name
        ))

    tiff_df = tiff_df %>% group_by(assay, tiff_file) %>% reframe(name = strsplit(name, " ")[[1]])
    tiff_df = tiff_df %>%
        mutate(name = ifelse(
            grepl("CTEBV[0-9]", name),
            sub("CTEBV", "CTEBV_", name),
            name
        ))
    .cp_name = function(x){
        x = sub("_2$", "", x)
        x = sub("Chiung", "Chung", x)
        cells = sapply(strsplit(x, "_"), function(x)x[length(x)])
        is_control = grepl("NegCTL", x) | grepl("Control", x)
        x[is_control]
        x[!is_control]
        control_group = ifelse(is_control, "NegCTL", "PosCTL")
        paste(cells, control_group, sep = '_')
    }
    tiff_df = tiff_df %>%
        mutate(name = ifelse(
            grepl("CellPellet", name),
            .cp_name(name),
            name
        ))

    tiff_df = tiff_df %>%
        mutate(name = sub("2468_", "", name)) %>%
        mutate(name = sub("_Rescanned_Cropped", "", name))

    tiff_df$sample_id = tiff_df$name
    tiff_df
}
