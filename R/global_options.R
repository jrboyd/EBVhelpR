#' @export
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
