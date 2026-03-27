# Internal IO helpers for phenocycler / RNAscope import workflows.

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
  all_dt_l <- list()
  for (name in names(files)) {
    message("Loading data for ", name)
    file_path <- files[[name]]
    stopifnot(file.exists(file_path))
    all_dt_l[[name]] <- .load_csv(file_path)
  }
  all_dt_l
}

#' Title
#'
#' @returns Character scalar path to the discovered original cell data directory.
#' @export
#'
#' @examples
#' \dontrun{
#' get_original_cell_data_dir()
#' }
get_original_cell_data_dir <- function() {
  env_dir <- Sys.getenv("EBVHELPER_DATA_DIR", unset = "")
  win_dir <- "C:/Users/boydj/OneDrive - UVM Larner College of Medicine/Lee, Kyra C's files - VolaricDataAndScriptsForJoe/"
  lin_dir <- "/gpfs1/home/j/r/jrboyd/VolaricDataAndScriptsForJoe/"

  candidates <- c(env_dir, win_dir, lin_dir)
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

.get_status_file <- function() {
  def_file <- "/gpfs1/pi/avolaric/files_jrboyd/EBVhelpR/inst/extdata/eber_status.xlsx"
  if (file.exists(def_file)) {
    return(def_file)
  }

  pkg_file <- system.file("extdata", "eber_status.xlsx", package = "EBVhelpR")
  if (!nzchar(pkg_file) || !file.exists(pkg_file)) {
    stop("Could not locate extdata/eber_status.xlsx.", call. = FALSE)
  }

  pkg_file
}

#' Title
#'
#' @returns Character scalar path to the package wrangled cell data directory.
#' @export
#'
#' @examples
#' \dontrun{
#' get_wrangled_cell_data_dir()
#' }
get_wrangled_cell_data_dir <- function() {
  file.path(get_original_cell_data_dir(), "../EBVhelpR_data")
}

