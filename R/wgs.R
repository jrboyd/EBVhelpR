.default_status_palette <- function(status) {
  status <- unique(as.character(stats::na.omit(status)))
  if (!length(status)) {
    return(character(0))
  }

  cols <- grDevices::hcl.colors(length(status), palette = "Dark 3")
  names(cols) <- status
  if ("need info" %in% names(cols)) {
    cols["need info"] <- "gray"
  }
  cols
}

.harmonize_bigwig_sample_ids <- function(sample, known_ids) {
  out <- as.character(sample)
  known_ids <- unique(as.character(known_ids))

  miss <- !(out %in% known_ids)
  out2 <- sub("_2$", "", out)
  out[miss & out2 %in% known_ids] <- out2[miss & out2 %in% known_ids]

  miss <- !(out %in% known_ids)
  out1 <- sub("_1$", "", out)
  out[miss & out1 %in% known_ids] <- out1[miss & out1 %in% known_ids]

  out
}

#' Load WGS reference ranges
#'
#' Reads the package chrSizes.txt file.
#'
#' [GenomicRanges::GRanges] object with sequence lengths set.
#'
#' @return A [GenomicRanges::GRanges] object.
#' @import Rsamtools
#' @import GenomicRanges
#' @export
load_wgs_reference_genome <- function() {
    chr_file = system.file("extdata/WGS_chrSizes.txt", package = "EBVhelpR", mustWork = TRUE)
    gr_df = read.table(chr_file)
    gr_df = gr_df %>% dplyr::rename(seqnames = V1, end = V2)
    gr_df$start = 1

  genome_gr <- GenomicRanges::GRanges(gr_df)
  sl <- GenomicRanges::width(genome_gr)
  names(sl) <- as.character(GenomicRanges::seqnames(genome_gr))
  sl <- sl[names(GenomeInfoDb::seqlengths(genome_gr))]
  GenomeInfoDb::seqlengths(genome_gr) <- sl

  genome_gr
}

#' Summarize host vs viral WGS read counts
#'
#' Loads per-sample count files and computes host/viral read fractions and
#' viral enrichment, matching the original WGS ranking workflow.
#'
#' @param count_dir Directory containing per-sample count txt files.
#' @param genome_gr Genome ranges from [load_wgs_reference_genome()].
#' @param viral_seqname Viral reference seqname in BAM/count files.
#' @param meta_df Optional metadata table with columns `sample_id` and
#'   `EBER_status`. Defaults to [load_meta_data()].
#' @param count_pattern Filename pattern for count files.
#'
#' @return A data frame ordered by increasing viral enrichment.
#' @export
load_wgs_count_summary <- function(
  count_dir,
  genome_gr,
  viral_seqname = "NC_007605.1",
  meta_df = NULL,
  count_pattern = "txt$"
) {
  count_files <- list.files(count_dir, pattern = count_pattern, full.names = TRUE)
  if (!length(count_files)) {
    stop("No count files found in `count_dir`.", call. = FALSE)
  }

  sample_id <- sub("_read.+", "", basename(count_files))
  count_df <- do.call(
    rbind,
    lapply(seq_along(count_files), function(i) {
      df <- utils::read.table(count_files[[i]], col.names = c("seqnames", "reads"))
      df$sample_id <- sample_id[[i]]
      df
    })
  )

  count_df <- dplyr::filter(count_df, .data$seqnames != "*")
  seq_lengths <- GenomeInfoDb::seqlengths(genome_gr)
  count_df$seqlengths <- as.numeric(seq_lengths[as.character(count_df$seqnames)])

  count_summary <- count_df |>
    dplyr::mutate(source = ifelse(.data$seqnames == viral_seqname, "viral", "host")) |>
    dplyr::group_by(.data$source, .data$sample_id) |>
    dplyr::summarise(
      total_reads = sum(.data$reads),
      total_length = sum(.data$seqlengths),
      .groups = "drop"
    ) |>
    tidyr::pivot_wider(
      id_cols = c("sample_id"),
      names_from = "source",
      values_from = c("total_reads", "total_length")
    ) |>
    dplyr::mutate(
      total_reads_viral = dplyr::coalesce(.data$total_reads_viral, 0),
      total_reads_host = dplyr::coalesce(.data$total_reads_host, 0),
      total_length_viral = dplyr::coalesce(.data$total_length_viral, 0),
      total_length_host = dplyr::coalesce(.data$total_length_host, 0),
      viral_length_fraction = .data$total_length_viral /
        (.data$total_length_host + .data$total_length_viral),
      viral_read_fraction = .data$total_reads_viral /
        (.data$total_reads_host + .data$total_reads_viral),
      viral_enrichment = .data$viral_read_fraction / .data$viral_length_fraction
    ) |>
    dplyr::arrange(.data$viral_enrichment)

  count_summary$sample_id <- factor(
    count_summary$sample_id,
    levels = unique(count_summary$sample_id)
  )

  if (is.null(meta_df)) {
    meta_df <- load_meta_data()
  }

  out <- merge(count_summary, meta_df, by = "sample_id", all.x = TRUE)
  out <- dplyr::mutate(
    out,
    EBER_status = ifelse(is.na(.data$EBER_status), "need info", .data$EBER_status)
  )
  out
}

#' Load and smooth EBV bigwig pileup profiles
#'
#' Fetches normalized bigwig signal across the EBV reference and merges status
#' metadata for downstream plotting.
#'
#' @param bigwig_dir Directory containing normalized bigwig files.
#' @param genome_gr Genome ranges from [load_wgs_reference_genome()].
#' @param viral_seqname Viral reference seqname to profile.
#' @param meta_df Optional metadata table with columns `sample_id` and
#'   `EBER_status`. Defaults to [load_meta_data()].
#' @param bigwig_pattern Pattern for normalized bigwig files.
#' @param smooth_n Moving-average window size.
#' @param mc_cores Number of cores to pass to seqsetvis.
#'
#' @return Data frame with columns including `sample`, `x`, `y`, and
#'   `EBER_status`.
#' @import seqsetvis
#' @export
load_wgs_bigwig_pileup <- function(
  bigwig_dir,
  genome_gr,
  viral_seqname = "NC_007605.1",
  meta_df = NULL,
  bigwig_pattern = "norm",
  smooth_n = 50,
  mc_cores = 1
) {
  if (is.null(meta_df)) {
    meta_df <- load_meta_data()
  }

  bw_files <- list.files(bigwig_dir, pattern = bigwig_pattern, full.names = TRUE)
  if (!length(bw_files)) {
    stop("No bigwig files found in `bigwig_dir`.", call. = FALSE)
  }

  bw_df <- data.frame(file = bw_files, stringsAsFactors = FALSE)
  bw_df$sample <- sub("_dedup.+", "", basename(bw_df$file))
  bw_df$sample <- .harmonize_bigwig_sample_ids(bw_df$sample, meta_df$sample_id)

  qgr <- genome_gr[as.character(GenomicRanges::seqnames(genome_gr)) == viral_seqname]
  if (!length(qgr)) {
    stop("`viral_seqname` was not found in `genome_gr`.", call. = FALSE)
  }

  old_cores <- getOption("mc.cores")
  on.exit(options(mc.cores = old_cores), add = TRUE)
  options(mc.cores = mc_cores)

  pileup_dt <- seqsetvis::ssvFetchBigwig(bw_df, qgr, return_data.table = TRUE)
  pileup_dt <- seqsetvis::applyMovingAverage(pileup_dt, n = smooth_n)

  bw_meta_df <- dplyr::filter(meta_df, .data$sample_id %in% bw_df$sample) |>
    dplyr::select(sample = .data$sample_id, .data$EBER_status)

  plot_dt <- merge(pileup_dt, bw_meta_df, by = "sample", all.x = TRUE)
  plot_dt$EBER_status[is.na(plot_dt$EBER_status)] <- "need info"
  plot_dt$x <- plot_dt$x - min(plot_dt$x)
  plot_dt
}

#' Plot WGS viral enrichment by sample
#'
#' @param wgs_count_summary Output from [load_wgs_count_summary()].
#' @param status_colors Named vector of colors keyed by `EBER_status`.
#'
#' @return A [ggplot2::ggplot] object.
#' @export
plot_wgs_viral_enrichment <- function(
  wgs_count_summary,
  status_colors = NULL
) {
  if (is.null(status_colors)) {
    status_colors <- .default_status_palette(wgs_count_summary$EBER_status)
  }

  ggplot2::ggplot(
    wgs_count_summary,
    ggplot2::aes(
      x = .data$sample_id,
      y = .data$viral_enrichment,
      color = .data$EBER_status,
      fill = .data$EBER_status
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::scale_color_manual(values = status_colors) +
    ggplot2::scale_fill_manual(values = status_colors) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "")
}

#' Plot WGS viral read totals by sample
#'
#' @param wgs_count_summary Output from [load_wgs_count_summary()].
#' @param status_colors Named vector of colors keyed by `EBER_status`.
#'
#' @return A [ggplot2::ggplot] object.
#' @export
plot_wgs_viral_reads <- function(
  wgs_count_summary,
  status_colors = NULL
) {
  if (is.null(status_colors)) {
    status_colors <- .default_status_palette(wgs_count_summary$EBER_status)
  }

  ggplot2::ggplot(
    wgs_count_summary,
    ggplot2::aes(
      x = .data$sample_id,
      y = .data$total_reads_viral,
      color = .data$EBER_status,
      fill = .data$EBER_status
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::scale_color_manual(values = status_colors) +
    ggplot2::scale_fill_manual(values = status_colors) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "")
}

#' Plot EBV bigwig pileup line profiles by status
#'
#' @param pileup_df Output from [load_wgs_bigwig_pileup()].
#' @param ylim Numeric length-2 vector giving y limits.
#'
#' @return A [ggplot2::ggplot] object.
#' @export
plot_wgs_pileup_lines <- function(
  pileup_df,
  ylim = c(0, 2000)
) {
  ggplot2::ggplot(
    pileup_df,
    ggplot2::aes(
      x = .data$x,
      y = .data$y,
      color = .data$sample,
      group = .data$sample
    )
  ) +
    ggplot2::geom_path() +
    ggplot2::facet_grid(.data$EBER_status ~ .) +
    ggplot2::guides(color = "none") +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::labs(y = "normalized read pileup", x = "position in EBV genome")
}

#' Plot EBV bigwig pileup heatmap
#'
#' @param pileup_df Output from [load_wgs_bigwig_pileup()].
#' @param wgs_count_summary Optional output from [load_wgs_count_summary()] used
#'   to order samples.
#' @param status_colors Named vector of colors keyed by `EBER_status`.
#' @param normalize_by_sample Logical; if `TRUE`, fill shows y / max(y) within
#'   each sample.
#' @param add_status_bar Logical; if `TRUE`, add a left-side annotation bar for
#'   status.
#'
#' @return A [ggplot2::ggplot] object.
#' @export
plot_wgs_pileup_heatmap <- function(
  pileup_df,
  wgs_count_summary = NULL,
  status_colors = NULL,
  normalize_by_sample = FALSE,
  add_status_bar = TRUE
) {
  plot_df <- pileup_df

  if (!is.null(wgs_count_summary)) {
    lvl <- levels(wgs_count_summary$sample_id)
    if (!is.null(lvl)) {
      plot_df$sample <- factor(plot_df$sample, levels = lvl)
    }
  }

  if (is.null(status_colors)) {
    status_colors <- .default_status_palette(plot_df$EBER_status)
  }

  if (normalize_by_sample) {
    fill_var <- "y_norm"
    plot_df <- plot_df |>
      dplyr::group_by(.data$sample) |>
      dplyr::mutate(y_norm = .data$y / max(.data$y, na.rm = TRUE)) |>
      dplyr::ungroup()
    fill_lab <- "sample normalized read pileup"
  } else {
    fill_var <- "y"
    fill_lab <- "depth normalized read pileup"
  }

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$x, y = .data$sample, fill = .data[[fill_var]])
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_x_continuous(labels = function(x) x / 1000) +
    ggplot2::labs(x = "EBV genomic position (kb)", y = "", fill = fill_lab)

  if (add_status_bar) {
    anno_dt <- unique(plot_df[, c("sample", "EBER_status")])
    anno_dt <- anno_dt[order(anno_dt$sample), , drop = FALSE]
    x_max <- max(plot_df$x)
    anno_dt$xmin <- -x_max * 0.07
    anno_dt$xmax <- -x_max * 0.025
    anno_dt$ymin <- as.numeric(anno_dt$sample) - 0.5
    anno_dt$ymax <- as.numeric(anno_dt$sample) + 0.5

    p <- p + ggplot2::annotate(
      "rect",
      xmin = anno_dt$xmin,
      xmax = anno_dt$xmax,
      ymin = anno_dt$ymin,
      ymax = anno_dt$ymax,
      fill = status_colors[anno_dt$EBER_status],
      color = status_colors[anno_dt$EBER_status]
    )
  }

  p
}
