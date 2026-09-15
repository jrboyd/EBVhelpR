#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import rlang
#' @importFrom glue glue
#' @importFrom lifecycle deprecated
#' @importFrom tibble tibble
## usethis namespace: end
NULL

# Columns referenced by non-standard evaluation inside dplyr/tidyr pipelines,
# and the magrittr `.` placeholder. Declared so R CMD check does not report them
# as undefined globals.
utils::globalVariables(c(
  ".",
  "Image Tag", "Sample", "SampleNumber", "V1", "XMax", "YMax",
  "assay", "name", "probe_control", "sample_id", "tiff_file"
))
