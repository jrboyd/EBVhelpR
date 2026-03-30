#' Supported Assay Type Constants
#'
#' Named list of supported assay type identifiers used by data loaders,
#' query constructors, and filtering helpers throughout the package.
#'
#' @format A named list with three character elements:
#' \describe{
#'   \item{phenocycler}{`"Phenocycler"`}
#'   \item{rnascope_4plex}{`"RNAScope_4plex"`}
#'   \item{rnascope_3plex+IF}{`"RNAScope_3plex+IF"`}
#' }
#' @examples
#' EBV_ASSAY_TYPES$rnascope_4plex
#' @export
EBV_ASSAY_TYPES = list(
    phenocycler = "Phenocycler",
    rnascope_4plex = "RNAScope_4plex",
    "rnascope_3plex+IF" = "RNAScope_3plex+IF"
)

assay_to_project_name = c(
    Phenocycler = "Phenocycler",
    RNAScope_4plex = "RNAScopeRound1",
    "RNAScope_3plex+IF" = "RNAScopeRound2"
)

project_name_to_assay = names(assay_to_project_name)
names(project_name_to_assay) = assay_to_project_name



.get_valid_project_names = function(){
    assay_to_project_name[unlist(EBV_ASSAY_TYPES)]
}
