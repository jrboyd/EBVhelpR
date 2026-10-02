# Data registry, resolver and local mirror.
#
# Every input the project reads is registered in inst/extdata/data_sources.csv,
# with its path relative to each storage location in
# inst/extdata/data_locations.csv. Both are plain CSV so they can be read and
# edited in Excel. Locations are ranked by `priority`: 1 is the most upstream
# copy (verified 2026-10-02 from creation times and name/size/mtime matching).
#
# A local mirror ($EBVHELPER_MIRROR_DIR) holds verbatim copies under
# cache/<mirror_subdir>/<relative path>, so the path itself records where a copy
# came from, and an append-only MANIFEST.tsv records size, mtime and md5.
#
# Nothing here moves or deletes data. Fetching only copies, and refuses to
# overwrite. Copies found to be redundant or corrupt are appended to a cleanup
# log for a person to act on.

.registry_csv <- function(name) {
  f <- system.file("extdata", name, package = "EBVhelpR", mustWork = TRUE)
  utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE,
                  na.strings = "", colClasses = "character")
}

#' Storage locations for project data, in priority order
#'
#' One row per location and path candidate. A location usually has several
#' candidates (UNC path, mapped drive, WSL mount, cluster path); the first that
#' exists on this machine is used. Candidates may contain `*` wildcards, which
#' are expanded with [Sys.glob()], so per-user paths need no editing.
#'
#' Read from `inst/extdata/data_locations.csv`, which is plain CSV so it can
#' be reviewed and edited in Excel.
#'
#' @returns A data frame with columns `location_id`, `priority` (integer, 1 =
#'   most upstream), `role` (`source`, `copy` or `derived`), `mirror_subdir`,
#'   `candidate` and `description`.
#' @seealso [ebv_data_sources()], [ebv_location_root()]
#' @export
#' @examples
#' unique(ebv_data_locations()[, c("location_id", "priority", "role")])
ebv_data_locations <- function() {
  x <- .registry_csv("data_locations.csv")
  x$priority <- as.integer(x$priority)
  x[order(x$priority), , drop = FALSE]
}

#' Registered project datasets
#'
#' One row per dataset. The columns named after a location (see
#' [ebv_data_locations()]) hold the dataset's path relative to that
#' location's root, blank where the location has no copy. `"."` means the
#' location root itself.
#'
#' `tier` says how a dataset is handled offline: `mirror` (copy it), `split`
#' (too big to keep whole; keep a per-sample reduction) or `reference` (never
#' copied by default; images and raw scans).
#'
#' @returns A data frame with columns `dataset_id`, `kind` (`file`/`dir`),
#'   `tier`, `assay`, `description`, one column per location, `bytes`,
#'   `mtime_utc` (of the highest-priority copy when registered) and `notes`.
#' @seealso [ebv_path()], [ebv_mirror_fetch()]
#' @export
#' @examples
#' head(ebv_data_sources()[, c("dataset_id", "tier", "labshare")])
ebv_data_sources <- function() {
  x <- .registry_csv("data_sources.csv")
  x$bytes <- as.numeric(x$bytes)
  locs <- unique(ebv_data_locations()$location_id)
  missing <- setdiff(locs, names(x))
  if (length(missing)) {
    stop("data_sources.csv lacks a column for location(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  x
}

.location_ids <- function() unique(ebv_data_locations()$location_id)

.dataset_row <- function(dataset_id) {
  src <- ebv_data_sources()
  row <- src[src$dataset_id == dataset_id, , drop = FALSE]
  if (nrow(row) != 1) {
    stop("Unknown dataset_id '", dataset_id, "'. See ebv_data_sources().", call. = FALSE)
  }
  row
}

#' Root directory of a storage location on this machine
#'
#' Checks `EBVHELPER_ROOT_<LOCATION_ID>` (upper case) first, then each
#' registered candidate in order, expanding `*` wildcards.
#'
#' @param location_id A `location_id` from [ebv_data_locations()].
#' @param must_work If `TRUE`, error when no candidate exists.
#' @returns The first existing directory, or `NA_character_`.
#' @export
#' @examples
#' ebv_location_root("labshare")
ebv_location_root <- function(location_id, must_work = FALSE) {
  locs <- ebv_data_locations()
  if (!location_id %in% locs$location_id) {
    stop("Unknown location_id '", location_id, "'.", call. = FALSE)
  }
  env <- Sys.getenv(paste0("EBVHELPER_ROOT_", toupper(location_id)), unset = "")
  cands <- c(if (nzchar(env)) env, locs$candidate[locs$location_id == location_id])
  cands <- unlist(lapply(cands, function(p) if (grepl("[*?]", p)) Sys.glob(p) else p))
  hit <- cands[dir.exists(cands)]
  if (length(hit)) return(hit[[1]])
  if (must_work) {
    stop("Location '", location_id, "' is not reachable from this machine. Tried:\n  ",
         paste(cands, collapse = "\n  "),
         "\nSet EBVHELPER_ROOT_", toupper(location_id), " to override.", call. = FALSE)
  }
  NA_character_
}

#' Local data mirror directory
#'
#' `EBVHELPER_MIRROR_DIR` if set, otherwise `~/EBV_data_mirror`. The
#' directory need not exist until something is fetched.
#'
#' @returns Character scalar path.
#' @export
ebv_mirror_dir <- function() {
  d <- Sys.getenv("EBVHELPER_MIRROR_DIR", unset = "")
  if (!nzchar(d)) d <- file.path(path.expand("~"), "EBV_data_mirror")
  d
}

.mirror_path <- function(location_id, rel) {
  sub <- unique(ebv_data_locations()$mirror_subdir[ebv_data_locations()$location_id == location_id])
  file.path(ebv_mirror_dir(), "cache", sub, rel)
}

.manifest_file <- function() file.path(ebv_mirror_dir(), "MANIFEST.tsv")

#' The mirror's fetch manifest
#'
#' Append-only record of every copy made into the mirror, by this package or by
#' `inst/scripts/mirror_fetch.ps1`.
#'
#' @returns A data frame (`source`, `local`, `bytes`, `source_mtime`, `md5`,
#'   `fetched`, `machine`), empty if nothing has been fetched.
#' @export
ebv_mirror_manifest <- function() {
  f <- .manifest_file()
  cols <- c("source", "local", "bytes", "source_mtime", "md5", "fetched", "machine")
  if (!file.exists(f)) {
    return(stats::setNames(data.frame(matrix(character(), 0, length(cols))), cols))
  }
  x <- utils::read.delim(f, stringsAsFactors = FALSE, colClasses = "character",
                         fileEncoding = "UTF-8-BOM", quote = "")
  x$bytes <- as.numeric(x$bytes)
  x
}

# Windows and WSL spellings of one mirror path must compare equal.
.norm_path <- function(p) {
  p <- gsub("\\\\", "/", p)
  p <- sub("^/mnt/([a-z])/", "\\U\\1:/", p, perl = TRUE)
  tolower(p)
}

.is_all_nul <- function(path, n = 512L) {
  if (!file.exists(path) || dir.exists(path) || file.size(path) == 0) return(FALSE)
  con <- file(path, "rb"); on.exit(close(con))
  b <- readBin(con, "raw", n = n)
  length(b) > 0 && all(b == as.raw(0))
}

# A mirror copy is trusted if it matches the size the manifest recorded for it
# (or, with no manifest row, the registry's size) and is not NUL-filled.
.FETCH_DONE <- ".ebv_fetch_complete"
.mirror_ok <- function(path, row) {
  # a directory counts only once a fetch of it has finished cleanly
  if (row$kind == "dir") return(file.exists(file.path(path, .FETCH_DONE)))
  if (!file.exists(path) || .is_all_nul(path)) return(FALSE)
  man <- ebv_mirror_manifest()
  rec <- man[.norm_path(man$local) == .norm_path(path), , drop = FALSE]
  want <- if (nrow(rec)) utils::tail(rec$bytes, 1) else row$bytes
  is.na(want) || file.size(path) == want
}

#' Candidate copies of a dataset on this machine
#'
#' Lists the mirror path and the source path for every location that registers
#' a copy, with whether each exists here. For inspection; [ebv_path()] makes
#' the choice.
#'
#' @param dataset_id A `dataset_id` from [ebv_data_sources()].
#' @returns A data frame with columns `location_id`, `priority`, `where`
#'   (`mirror`/`source`), `path`, `exists`.
#' @export
ebv_dataset_paths <- function(dataset_id) {
  row  <- .dataset_row(dataset_id)
  locs <- unique(ebv_data_locations()[, c("location_id", "priority")])
  out <- lapply(seq_len(nrow(locs)), function(i) {
    id <- locs$location_id[i]; rel <- row[[id]]
    if (is.na(rel)) return(NULL)
    root <- ebv_location_root(id)
    data.frame(location_id = id, priority = locs$priority[i],
               where = c("mirror", "source"),
               path = c(.mirror_path(id, rel), if (is.na(root)) NA_character_ else file.path(root, rel)),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  out$exists <- !is.na(out$path) & (file.exists(out$path) | dir.exists(out$path))
  out
}

#' Resolve a registered dataset to a path on this machine
#'
#' Order: a verified mirror copy (any location, highest priority first), then
#' the highest-priority reachable source. Fails with a message naming every
#' path tried when neither is available, so a script never silently reads a
#' stale or partial copy.
#'
#' @param dataset_id A `dataset_id` from [ebv_data_sources()].
#' @param use_mirror Consider mirror copies (default `TRUE`).
#' @param must_work If `FALSE`, return `NA_character_` instead of erroring.
#' @returns Character scalar path.
#' @seealso [ebv_mirror_fetch()], [ebv_dataset_paths()]
#' @export
#' @examples
#' \dontrun{
#' ihc <- read.csv(ebv_path("ihc_clean_2025_09_03"))
#' }
ebv_path <- function(dataset_id, use_mirror = TRUE, must_work = TRUE) {
  row   <- .dataset_row(dataset_id)
  cands <- ebv_dataset_paths(dataset_id)
  if (use_mirror) {
    for (p in cands$path[cands$where == "mirror" & cands$exists]) {
      if (.mirror_ok(p, row)) return(p)
    }
  }
  src <- cands[cands$where == "source" & cands$exists, , drop = FALSE]
  if (nrow(src)) {
    p <- src$path[[1]]
    if (row$kind == "file" && .is_all_nul(p)) {
      warning("Source copy of '", dataset_id, "' at ", p, " is NUL-filled; trying lower-priority copies.",
              call. = FALSE)
      ok <- !vapply(src$path, .is_all_nul, logical(1))
      if (!any(ok)) p <- NA_character_ else p <- src$path[ok][[1]]
    }
    if (!is.na(p)) return(p)
  }
  if (!must_work) return(NA_character_)
  stop("Dataset '", dataset_id, "' is not available on this machine.\nTried:\n  ",
       paste(sprintf("[%s %s] %s", cands$location_id, cands$where,
                     ifelse(is.na(cands$path), "(location not reachable)", cands$path)), collapse = "\n  "),
       "\nWhile a source is reachable, run ebv_mirror_fetch(\"", dataset_id, "\").", call. = FALSE)
}

.append_manifest <- function(source, local, bytes, source_mtime, md5) {
  f <- .manifest_file()
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(f)) {
    cat("source\tlocal\tbytes\tsource_mtime\tmd5\tfetched\tmachine\n", file = f)
  }
  stamp <- function(t) format(t, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  cat(paste(source, local, format(bytes, scientific = FALSE), stamp(source_mtime), md5,
            stamp(Sys.time()), Sys.info()[["nodename"]], sep = "\t"), "\n",
      file = f, append = TRUE, sep = "")
}

.copy_one <- function(from, to, hash) {
  if (file.exists(to)) {
    if (file.size(to) == file.size(from) &&
        abs(as.numeric(difftime(file.mtime(to), file.mtime(from), units = "secs"))) <= 2) {
      return("unchanged")
    }
    ebv_log_cleanup(to, "mirror copy differs from its source; not overwritten (copy-only policy)",
                    keep = from, location = "mirror")
    return("differs")
  }
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(from, to, overwrite = FALSE, copy.date = TRUE)) stop("copy failed: ", from, call. = FALSE)
  if (file.size(to) != file.size(from)) stop("size mismatch after copying ", from, call. = FALSE)
  md5 <- if (hash) unname(tools::md5sum(to)) else ""
  .append_manifest(from, to, file.size(from), file.mtime(from), md5)
  "copied"
}

#' Copy a registered dataset into the local mirror
#'
#' Copy only: nothing at the source is touched, and an existing mirror copy is
#' never overwritten. A mirror copy that differs from its source is left in
#' place and logged with [ebv_log_cleanup()].
#'
#' For files of tens of GB on Windows, `inst/scripts/mirror_fetch.ps1`
#' (robocopy) is much faster than R reading through a WSL mount, and writes
#' the same manifest and mirror layout.
#'
#' @param dataset_id A `dataset_id` from [ebv_data_sources()].
#' @param from Location to copy from; default the highest-priority reachable one.
#' @param hash Record an md5 of each copied file (default `TRUE`).
#' @param allow_reference Datasets of tier `reference` (images, raw scans) are
#'   refused unless this is `TRUE`.
#' @returns Invisibly, a data frame of files and what happened to each
#'   (`copied`, `unchanged`, `differs`).
#' @export
ebv_mirror_fetch <- function(dataset_id, from = NULL, hash = TRUE, allow_reference = FALSE) {
  row <- .dataset_row(dataset_id)
  if (row$tier == "reference" && !allow_reference) {
    stop("'", dataset_id, "' is tier 'reference' (", row$description, "); not mirrored by default. ",
         "Use allow_reference = TRUE to copy it anyway.", call. = FALSE)
  }
  if (!is.na(row$notes) && grepl("Never mirror", row$notes, fixed = TRUE)) {
    stop("'", dataset_id, "' is marked 'Never mirror': ", row$notes, call. = FALSE)
  }
  cands <- ebv_dataset_paths(dataset_id)
  src <- cands[cands$where == "source" & cands$exists, , drop = FALSE]
  if (!is.null(from)) src <- src[src$location_id == from, , drop = FALSE]
  if (!nrow(src)) {
    stop("No reachable source for '", dataset_id, "'", if (!is.null(from)) paste0(" at '", from, "'"), ".",
         call. = FALSE)
  }
  loc <- src$location_id[[1]]; s <- src$path[[1]]
  d <- .mirror_path(loc, row[[loc]])
  pairs <- if (row$kind == "dir") {
    rel <- list.files(s, recursive = TRUE, all.files = TRUE, no.. = TRUE)
    rel <- rel[basename(rel) != .FETCH_DONE]
    data.frame(from = file.path(s, rel), to = file.path(d, rel), stringsAsFactors = FALSE)
  } else data.frame(from = s, to = d, stringsAsFactors = FALSE)
  pairs$result <- vapply(seq_len(nrow(pairs)), function(i) .copy_one(pairs$from[i], pairs$to[i], hash), "")
  if (row$kind == "dir" && all(pairs$result %in% c("copied", "unchanged"))) {
    writeLines(c(paste("source:", s), paste("files:", nrow(pairs)),
                 paste("completed:", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))),
               file.path(d, .FETCH_DONE))
  }
  message(sprintf("%s from %s: %d copied, %d unchanged, %d differ (logged)", dataset_id, loc,
                  sum(pairs$result == "copied"), sum(pairs$result == "unchanged"), sum(pairs$result == "differs")))
  invisible(pairs)
}

#' Cleanup log location
#'
#' `EBVHELPER_CLEANUP_LOG` if set, otherwise `cleanup_targets.tsv` in
#' [ebv_mirror_dir()].
#'
#' @returns Character scalar path.
#' @export
ebv_cleanup_log <- function() {
  f <- Sys.getenv("EBVHELPER_CLEANUP_LOG", unset = "")
  if (!nzchar(f)) f <- file.path(ebv_mirror_dir(), "cleanup_targets.tsv")
  f
}

#' Record a file or directory as a cleanup candidate
#'
#' Appends a row to the cleanup log; the target itself is never touched. A
#' repeat of the same `path` and `reason` is ignored, so audits can be re-run.
#'
#' @param path The redundant or corrupt copy.
#' @param reason Why it is a candidate (e.g. "truncated copy", "duplicate").
#' @param keep The copy to keep instead, if any.
#' @param location Where `path` lives (a `location_id`, `"mirror"`, ...).
#' @param bytes Size, if known (looked up when `path` is reachable).
#' @param log_file Defaults to [ebv_cleanup_log()].
#' @returns Invisibly, `TRUE` if a row was added.
#' @export
ebv_log_cleanup <- function(path, reason, keep = NA_character_, location = NA_character_,
                            bytes = NA_real_, log_file = ebv_cleanup_log()) {
  if (is.na(bytes) && file.exists(path) && !dir.exists(path)) bytes <- file.size(path)
  if (file.exists(log_file)) {
    old <- utils::read.delim(log_file, stringsAsFactors = FALSE, colClasses = "character", quote = "")
    if (any(old$path == path & old$reason == reason)) return(invisible(FALSE))
  } else {
    dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
    cat("logged_utc\tmachine\tlocation\tpath\tbytes\treason\tkeep\n", file = log_file)
  }
  clean <- function(x) gsub("[\t\r\n]", " ", ifelse(is.na(x), "", as.character(x)))
  cat(paste(format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), Sys.info()[["nodename"]],
            clean(location), clean(path), clean(if (is.na(bytes)) "" else format(bytes, scientific = FALSE)),
            clean(reason), clean(keep), sep = "\t"), "\n", file = log_file, append = TRUE, sep = "")
  invisible(TRUE)
}
