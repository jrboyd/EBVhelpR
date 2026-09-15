#' Supported Assay Type Constants
#'
#' Named list of supported assay type identifiers used by data loaders,
#' query constructors, and filtering helpers throughout the package.
#'
#' @format A named list of three character scalars, the assay type identifiers:
#' \describe{
#'   \item{Phenocycler}{`"Phenocycler"`}
#'   \item{RNAScope_4plex}{`"RNAScope_4plex"`}
#'   \item{RNAScope_3plex+IF}{`"RNAScope_3plex+IF"`}
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


#' Colors used for each EBER status
#'
#' Reads the packaged palette so that EBER status colors stay consistent across
#' every figure.
#'
#' @returns A named character vector of colors, named by EBER status.
#' @export
#'
#' @examples
#' get_colors_EBER_status()
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
    "Autofluorescence"
)
stopifnot(EBV_CHANNELS$RNAScope_4plex[6] == "Autofluorescence")

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
    "EBNA1-Ab",
    "Autofluorescence"
)
stopifnot(EBV_CHANNELS$`RNAScope_3plex+IF`[6] == "Autofluorescence")

#' Channel identities in tiff files.
#'
#' Named list of channel names, one character vector per assay. The vectors are
#' positional: element `n` is the identity of channel `n` in the TIFF, which is
#' how they are passed to `channel_names` when fetching crops.
#'
#' The autofluorescence channel of both RNAscope panels is named
#' `"Autofluorescence"`. Before version 0.1.8 it was misspelled
#' `"Autofluoresence"`, so code matching that string will silently stop
#' matching; the corrected spelling agrees with the analysis scripts.
#'
#' @format A named list of three character vectors, named by assay type:
#' \describe{
#'   \item{Phenocycler}{16 channels, DAPI first.}
#'   \item{RNAScope_4plex}{6 channels: DAPI, EBER1, EBNA2, LMP1, EBNA3, Autofluorescence.}
#'   \item{RNAScope_3plex+IF}{6 channels: DAPI, EBER, LMP1, EBNA1, EBNA1-Ab, Autofluorescence.}
#' }
#' @examples
#' EBV_CHANNELS$RNAScope_4plex
#' @export
EBV_CHANNELS = EBV_CHANNELS


EBV_OPAL_DECODE = EBV_ASSAY_TYPES
EBV_OPAL_DECODE$Phenocycler = NULL
EBV_OPAL_DECODE$RNAScope_4plex = c(
    Opal520 = "EBER",
    Opal620 = "EBNA2",
    Opal570 = "LMP1",
    Opal690 = "EBNA3"
)
EBV_OPAL_DECODE$`RNAScope_3plex+IF` = c(
    Opal520 = "EBER",
    Opal620 = "EBNA1",
    Opal570 = "LMP1",
    Opal690 = "EBNA1-Ab"
)


#' Opal dye to target decode for the RNAscope panels
#'
#' Named list, one character vector per RNAscope assay, whose names are Opal dye
#' strings and whose values are the targets those dyes label. `Phenocycler` is
#' absent, so `EBV_OPAL_DECODE[["Phenocycler"]]` is `NULL`.
#'
#' @section Unusable channels:
#' Not every dye listed here yields usable data. `Opal690` is present in the
#' images but is **not** in either per-cell export, so only three of the four
#' probes reach the user per cell. Check [EBV_OPAL_USABLE] before assuming a
#' target exists; the entries are kept here because the channels really are in
#' the images.
#'
#' The two failures differ. In `RNAScope_3plex+IF` (the EBNA1 antibody), Opal690
#' is the channel most correlated with the dedicated autofluorescence channel in
#' all three images tested (Spearman rho 0.33-0.52, partial on DAPI; two to ten
#' times any other Opal) and has no zero-pixel population, sitting on a floor of
#' 2.5-2.9 where Opal520/570 sit at exactly 0.00 — it is measuring
#' autofluorescence. In `RNAScope_4plex` (EBNA3) it is simply empty: dynamic
#' range 1.2x above background and rho -0.002 against autofluorescence, while
#' Opal520 (EBER1) in the same image gives 42x.
#'
#' @format A named list of two character vectors, each of length four, named by
#'   Opal dye:
#' \describe{
#'   \item{RNAScope_4plex}{`Opal520` EBER, `Opal620` EBNA2, `Opal570` LMP1,
#'     `Opal690` EBNA3 (unusable).}
#'   \item{RNAScope_3plex+IF}{`Opal520` EBER, `Opal620` EBNA1, `Opal570` LMP1,
#'     `Opal690` EBNA1-Ab (unusable).}
#' }
#' @seealso [EBV_OPAL_USABLE]
#' @examples
#' EBV_OPAL_DECODE$RNAScope_4plex
#' @export
EBV_OPAL_DECODE = EBV_OPAL_DECODE

EBV_OPAL_USABLE = lapply(EBV_OPAL_DECODE, function(x){
    stats::setNames(names(x) != "Opal690", names(x))
})

#' Which Opal channels carry usable signal
#'
#' Mirrors [EBV_OPAL_DECODE] in shape, with a logical per Opal dye. `Opal690` is
#' `FALSE` in both RNAscope panels — in `RNAScope_3plex+IF` it measures
#' autofluorescence, in `RNAScope_4plex` it is empty. See the Unusable channels
#' section of [EBV_OPAL_DECODE] for the supporting measurements.
#'
#' Use this rather than assuming every target in `EBV_OPAL_DECODE` reached the
#' per-cell export: `Opal690` did not.
#'
#' @format A named list of two logical vectors, named by Opal dye, matching the
#'   layout of [EBV_OPAL_DECODE].
#' @seealso [EBV_OPAL_DECODE]
#' @examples
#' # targets that actually carry signal for the 4plex panel
#' EBV_OPAL_DECODE$RNAScope_4plex[EBV_OPAL_USABLE$RNAScope_4plex]
#' @export
EBV_OPAL_USABLE = EBV_OPAL_USABLE
