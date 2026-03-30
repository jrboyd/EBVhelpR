.expand_previous_samples <- function(meta_df) {
  prev_sample <- grepl("previously", meta_df$sample_id)
  if (!any(prev_sample)) {
    return(meta_df)
  }

  rows <- meta_df[prev_sample, , drop = FALSE]
  expanded <- lapply(seq_len(nrow(rows)), function(i) {
    parts <- strsplit(rows$sample_id[i], " previously ", fixed = TRUE)[[1]]
    data.frame(sample_id = parts, EBER_status = rows$EBER_status[i], stringsAsFactors = FALSE)
  })

  rbind(meta_df[!prev_sample, , drop = FALSE], do.call(rbind, expanded))
}

.expand_slash_samples <- function(meta_df) {
  split_ids <- grepl("/", meta_df$sample_id)
  if (!any(split_ids)) {
    return(meta_df)
  }

  split_df <- meta_df[split_ids, , drop = FALSE]
  sp <- strsplit(split_df$sample_id, "/", fixed = TRUE)
  left_df <- split_df
  right_df <- split_df
  left_df$sample_id <- vapply(sp, function(x) x[1], character(1))
  right_df$sample_id <- vapply(sp, function(x) x[2], character(1))

  rbind(meta_df[!split_ids, , drop = FALSE], left_df, right_df)
}

#' Load full EBER clinical metadata
#'
#' Reads the packaged EBER annotation file, standardizes sample_id naming, and
#' expands composite entries into one row per sample_id.
#'
#' @return A data frame with columns `sample_id` and `EBER_status`.
#' @examples
#' \dontrun{
#' meta_df <- load_meta_data()
#' }
#' @export
load_meta_data <- function() {
  meta_df <- openxlsx::read.xlsx(.get_status_file())

  if (ncol(meta_df) < 2) {
    stop("EBER status sheet must have at least two columns.", call. = FALSE)
  }

  meta_df <- meta_df[, 1:2, drop = FALSE]
  colnames(meta_df) <- c("sample_id", "EBER_status")

  meta_df$sample_id <- gsub("-", "_", meta_df$sample_id)
  meta_df$sample_id <- gsub("[()]", "", meta_df$sample_id)

  meta_df <- .expand_previous_samples(meta_df)
  meta_df <- .expand_slash_samples(meta_df)

  meta_df
}