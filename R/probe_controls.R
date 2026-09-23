#' Probe-control sections and the sample ids they share
#'
#' Some sections were stained with the vendor control cocktail instead of the
#' EBV panel: RNAscope `*_PosCTL` / `*_NegCTL` serial sections and
#' `301_CellPellet_NegCTL_*` pellets, and the Phenocycler
#' `CellPelletSlide_Control_*` slide. The loaders flag these in
#' `probe_control` and give them the SAME `sample_id` as the sample they were
#' cut from, so that a control can be linked back to its sample. The cost is
#' that anything aggregating by `sample_id` silently pools a control section
#' into its sample. In P1 that put CTEBV_15's positive-probe section into the
#' Latency II anchor, where it supplied 97% of the "EBER1+" cells.
#'
#' Use `unique_id` (which appends `probe_control`), or filter
#' `probe_control == ""`, before grouping by sample. The loaders warn when their
#' output contains probe-control rows; set
#' `options(EBVhelpR.warn_probe_controls = FALSE)` to silence that.
#'
#' @name probe_controls
NULL

## Matches the suffixes the collaborator exports use. Deliberately NOT a bare
## "[Pp]os" / "[Nn]eg", which would also flag any sample name containing them.
.PROBE_NEG_PATTERN <- "NegCTL"
.PROBE_POS_PATTERN <- "PosCTL"

.warn_probe_controls <- function(df, caller) {
    if (!isTRUE(getOption("EBVhelpR.warn_probe_controls", TRUE))) {
        return(invisible(df))
    }
    if (!nrow(df) || !"probe_control" %in% colnames(df)) {
        return(invisible(df))
    }
    is_probe <- !is.na(df$probe_control) & df$probe_control != ""
    if (!any(is_probe)) {
        return(invisible(df))
    }
    ids <- sort(unique(df$sample_id[is_probe]))
    warning(
        caller, "(): ", sum(is_probe), " row(s) are probe-control sections that share ",
        "sample_id with their parent sample (", paste(utils::head(ids, 8), collapse = ", "),
        if (length(ids) > 8) sprintf(", ... %d total", length(ids)) else "", "). ",
        "Grouping by sample_id pools them into the sample; group by unique_id or ",
        "filter probe_control == \"\" first. ",
        "Silence with options(EBVhelpR.warn_probe_controls = FALSE).",
        call. = FALSE
    )
    invisible(df)
}

.without_probe_warnings <- function(expr) {
    old <- options(EBVhelpR.warn_probe_controls = FALSE)
    on.exit(options(old), add = TRUE)
    force(expr)
}

#' Check RNAscope summaries against the per-cell files they should describe
#'
#' The collaborator summaries are keyed by a sample name, and a name can cover
#' more than the one section a caller expects. The 4plex cell-pellet summary
#' (`4PlexRNAScopeCellPelletSingleExpression_2026-01-12.csv`) merges each cell
#' line's negative-probe pellet into it and leaves no flag to filter on. D_EB_3's
#' 4plex summary sums both of its scans. Neither can be seen from the summary
#' alone. Both show up as a cell total that disagrees with the per-cell files.
#'
#' Totals are compared per `(assay, unique_id)`. `cell_count` in
#' [load_cell_source_files()] is a `wc -l` line count, so one header line per
#' file is subtracted.
#'
#' @param summary_df RNAscope summary rows, e.g. from
#'   [load_rnascope_summary_files()]. Needs `assay`, `sample_id`,
#'   `probe_control` and `count`.
#' @param cell_files_df Per-cell file table, e.g. from
#'   [load_cell_source_files()] or [get_query_cell_files_df()]. Needs `assay`,
#'   `sample_id`, `probe_control` and `cell_count`. Restrict it to the scans you
#'   intend to use: a summary covering a scan you dropped then reports as a
#'   mismatch, which is the point.
#' @param warn Warn when any total disagrees. Defaults to `TRUE`.
#'
#' @return A data frame with one row per `(assay, unique_id)`: `n_summary`,
#'   `n_cell_files`, and `match` (`NA` when present on only one side).
#' @examples
#' \dontrun{
#' chk <- check_summary_against_cell_files(
#'   load_rnascope_summary_files(),
#'   load_cell_source_files()
#' )
#' subset(chk, !match)
#' }
#' @export
check_summary_against_cell_files <- function(summary_df, cell_files_df, warn = TRUE) {
    need_s <- c("assay", "sample_id", "probe_control", "count")
    need_c <- c("assay", "sample_id", "probe_control", "cell_count")
    miss_s <- setdiff(need_s, colnames(summary_df))
    miss_c <- setdiff(need_c, colnames(cell_files_df))
    if (length(miss_s)) stop("summary_df lacks column(s): ", paste(miss_s, collapse = ", "), call. = FALSE)
    if (length(miss_c)) stop("cell_files_df lacks column(s): ", paste(miss_c, collapse = ", "), call. = FALSE)

    key <- function(df) ifelse(is.na(df$probe_control) | df$probe_control == "",
                               df$sample_id, paste(df$sample_id, df$probe_control))
    s <- data.frame(assay = summary_df$assay, unique_id = key(summary_df),
                    n = as.numeric(summary_df$count), stringsAsFactors = FALSE)
    s <- s[!is.na(s$n), , drop = FALSE]
    s <- stats::aggregate(n ~ assay + unique_id, data = s, FUN = sum)
    names(s)[names(s) == "n"] <- "n_summary"

    cf <- cell_files_df[cell_files_df$assay %in% unique(summary_df$assay), , drop = FALSE]
    c_ <- data.frame(assay = cf$assay, unique_id = key(cf),
                     n = as.numeric(cf$cell_count) - 1, stringsAsFactors = FALSE)
    c_ <- stats::aggregate(n ~ assay + unique_id, data = c_, FUN = sum)
    names(c_)[names(c_) == "n"] <- "n_cell_files"

    out <- merge(s, c_, by = c("assay", "unique_id"), all = TRUE)
    out$match <- out$n_summary == out$n_cell_files
    out <- out[order(out$assay, out$unique_id), , drop = FALSE]
    rownames(out) <- NULL

    bad <- out[!is.na(out$match) & !out$match, , drop = FALSE]
    if (warn && nrow(bad)) {
        warning(
            nrow(bad), " summary total(s) disagree with their per-cell files: ",
            paste(sprintf("%s:%s (%s vs %s)", bad$assay, bad$unique_id,
                          format(bad$n_summary, big.mark = ","),
                          format(bad$n_cell_files, big.mark = ",")), collapse = "; "),
            ". The summary covers a different set of sections; recount from the per-cell files.",
            call. = FALSE
        )
    }
    out
}
