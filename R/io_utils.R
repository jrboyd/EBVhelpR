# Internal IO helpers for Phenocycler / RNAscope import workflows.

.load_csv <- function(file_path) {
  tryCatch(
    {
      readr::read_csv(file_path, show_col_types = FALSE)
    },
    error = function(e) {
      message("Error loading file ", file_path, ": ", e$message)
      NULL
    }
  )
}

.load_csv_list <- function(files) {
  all_dt_l = lapply(files, function(file_path){
      message("Loading data for ", file_path)
      stopifnot(file.exists(file_path))
      out = .load_csv(file_path)
      out$source_file = file_path
      out
  })
  all_dt_l
}

#' Get Original Cell Data Directory
#'
#' Resolves the root directory containing raw cell data inputs by checking
#' `EBVHELPER_DATA_DIR` first, then the `kyra_onedrive` location of the data
#' registry ([ebv_location_root()]), then the legacy default paths.
#'
#' @return Character scalar path to the discovered original cell data directory.
#' @examples
#' \dontrun{
#' get_original_cell_data_dir()
#' }
#' @export
get_original_cell_data_dir <- function() {
  env_dir <- Sys.getenv("EBVHELPER_DATA_DIR", unset = "")
  reg_dir <- ebv_location_root("kyra_onedrive")
  win_dir <- "C:/Users/boydj/OneDrive - UVM Larner College of Medicine/Lee, Kyra C's files - VolaricDataAndScriptsForJoe/"
  lin_dir <- "/gpfs1/home/j/r/jrboyd/VolaricDataAndScriptsForJoe/"

  candidates <- c(env_dir, if (!is.na(reg_dir)) reg_dir, win_dir, lin_dir)
  candidates <- candidates[nzchar(candidates)]
  existing <- candidates[dir.exists(candidates)]

  if (!length(existing)) {
    stop(
      "No data directory found. Set EBVHELPER_DATA_DIR or ensure one of the default paths exists.",
      call. = FALSE
    )
  }

  existing[[1]]
}

#' Locate the EBER status sheet
#'
#' Prefers the copy installed with the package so that the same code gives the
#' same EBER calls everywhere. An alternative sheet can be selected explicitly
#' with `EBVHELPER_STATUS_FILE`; this used to happen implicitly via a hardcoded
#' cluster path, which silently shadowed the packaged copy on the VACC.
#'
#' @return Character scalar path to the EBER status workbook.
#' @noRd
.get_status_file <- function() {
  env_file <- Sys.getenv("EBVHELPER_STATUS_FILE", unset = "")
  if (nzchar(env_file)) {
    if (!file.exists(env_file)) {
      stop(
        "EBVHELPER_STATUS_FILE is set to a file that does not exist: ", env_file,
        call. = FALSE
      )
    }
    message("Using EBER status file from EBVHELPER_STATUS_FILE: ", env_file)
    return(env_file)
  }

  pkg_file <- system.file("extdata", "eber_status.xlsx", package = "EBVhelpR")
  if (!nzchar(pkg_file) || !file.exists(pkg_file)) {
    stop("Could not locate extdata/eber_status.xlsx.", call. = FALSE)
  }

  pkg_file
}

#' Get Wrangled Cell Data Directory
#'
#' Returns the directory where processed package-ready cell data files are
#' written. If `EBVHELPER_DATA_DIR` is set this is its sibling `EBVhelpR_data`,
#' as before. Otherwise it is the registered `cell_stores` dataset
#' ([ebv_path()]: a verified mirror copy, else the `ebvhelpr_data` location),
#' falling back to the sibling of [get_original_cell_data_dir()].
#'
#' @return Character scalar path to the package wrangled cell data directory.
#' @examples
#' \dontrun{
#' get_wrangled_cell_data_dir()
#' }
#' @export
get_wrangled_cell_data_dir <- function() {
  if (!nzchar(Sys.getenv("EBVHELPER_DATA_DIR", unset = ""))) {
    reg <- ebv_path("cell_stores", must_work = FALSE)
    if (!is.na(reg)) return(reg)
  }
  file.path(get_original_cell_data_dir(), "../EBVhelpR_data")
}

