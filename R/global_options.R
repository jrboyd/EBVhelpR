#' Supported Assay Type Constants
#'
#' Named list of supported assay type identifiers used by data loaders,
#' query constructors, and filtering helpers throughout the package.
#'
#' @format A named list with three character elements:
#' \describe{
#'   \item{Phenocycler}{`"Phenocycler"`}
#'   \item{RNAScope_4plex}{`"RNAScope_4plex"`}
#'   \item{rnascope_3plex+IF}{`"RNAScope_3plex+IF"`}
#' }
#' @examples
#' EBV_ASSAY_TYPES$RNAScope_4plex
#' @export
EBV_ASSAY_TYPES = list(
    Phenocycler = "Phenocycler",
    RNAScope_4plex = "RNAScope_4plex",
    "RNAScope_3plex+IF" = "RNAScope_3plex+IF"
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


#' Title
#'
#' @returns
#' @export
#'
#' @examples
get_colors_EBER_status = function(){
    readRDS(system.file("extdata/colors_EBER_status.Rds", package = "EBVhelpR", mustWork = TRUE))
}


EBV_CHANNELS = EBV_ASSAY_TYPES

# 1
# DAPI
# 2
# PAX5
# 3
# EBNA2
# 4
# CD3
# 5
# EBNA3A
# 6
# LMP1
# 7
# CD4
# 8
# CD30
# 9
# c-Myc
# 10
# EBNA3B
# 11
# CD8
# 12
# EBNA3C
# 13
# PDL1
# 14
# LMP2A
# 15
# CD20
# 16
# PD1

EBV_CHANNELS$Phenocycler = c(
    "DAPI",
    'PAX5',
    'EBNA2',
    "CD3",
    "EBNA3A",
    "LMP1",
    "CD4",
    "CD30",
    "c-Myc",
    "EBNA3B",
    "CD8",
    "EBNA3C",
    "PDL1",
    "LMP2A",
    "CD20",
    "PD1"
)
stopifnot(EBV_CHANNELS$Phenocycler[6] == "LMP1")
stopifnot(EBV_CHANNELS$Phenocycler[11] == "CD8")
stopifnot(EBV_CHANNELS$Phenocycler[15] == "CD20")
stopifnot(EBV_CHANNELS$Phenocycler[16] == "PD1")


# Probe Cocktail
# C1 Dye
# C2 Dye
# C3 Dye
# C4 Dye
# Test Probe Cocktail (C1- EBER1/C2-EBNA2/C3-LMP1/C4-EBNA3)
# 520
# 620
# 570
# 690
# The other two channels are DAPI and autofluorescence.
#


EBV_CHANNELS$RNAScope_4plex = c(
    "DAPI",
    "EBER1",
    "EBNA2",
    "LMP1",
    "EBNA3",
    "Autofluoresence"
)

# For the RNAScope + IF the channels are set up like this:
#
#     
# Probe Cocktail
# C1 Dye (TSA F1)
# C3 Dye (TSA F2)
# C4 Dye (TSA F3)
# EBNA 1 Ab Dye (TSA F4)
# C1 -EBER/C3- LMP1/C4-EBNA1
# 520
# 570
# 620
# 690

EBV_CHANNELS$`RNAScope_3plex+IF` = c(
    "DAPI",
    "EBER",
    "LMP1",
    "EBNA1",
    "Autofluoresence"
)

#' Channel identities in tiff files.
#'
#' Named list of channel names.
#'
#' @format A named list with three character elements:
#' \describe{
#'   \item{Phenocycler}{`"Phenocycler"`}
#'   \item{RNAScope_4plex}{`"RNAScope_4plex"`}
#'   \item{rnascope_3plex+IF}{`"RNAScope_3plex+IF"`}
#' }
#' @examples
#' EBV_CHANNELS$RNAScope_4plex
#' @export
EBV_CHANNELS = EBV_CHANNELS
