.stop_if_missing_tiffplotr <- function() {
    if (!requireNamespace("TiffPlotR", quietly = TRUE)) {
        stop(
            "Package `TiffPlotR` is required for TIFF image extraction. Install it first.",
            call. = FALSE
        )
    }
}

.find_tiff_file_by_sample <- function(object, sample_id) {
    tiff_df <- get_query_tiff_paths_df(object)
    sel <- tiff_df[tiff_df$sample_id == sample_id, , drop = FALSE]

    if (!nrow(sel)) {
        stop("No matching TIFF found for sample_id: ", sample_id, call. = FALSE)
    }

    if (nrow(sel) > 1) {
        warning(
            "Multiple TIFF files found for sample_id ",
            sample_id,
            "; using the first one.",
            call. = FALSE
        )
    }

    sel$tiff_file[[1]]
}

.rects_from_df <- function(df) {
    lapply(seq_len(nrow(df)), function(i) {
        TiffPlotR::TiffRect(df$XMin[i], df$XMax[i], df$YMin[i], df$YMax[i])
    })
  if (!requireNamespace("TiffPlotR", quietly = TRUE)) {
    stop(
      "Package `TiffPlotR` is required for TIFF image extraction. Install it first.",
      call. = FALSE
    )
  }
}

.find_tiff_file_by_sample <- function(object, sample_id) {
  tiff_df <- get_query_tiff_paths_df(object)
  sel <- tiff_df[tiff_df$sample_id == sample_id, , drop = FALSE]

  if (!nrow(sel)) {
    stop("No matching TIFF found for sample_id: ", sample_id, call. = FALSE)
  }

  if (nrow(sel) > 1) {
    warning(
      "Multiple TIFF files found for sample_id ",
      sample_id,
      "; using the first one.",
      call. = FALSE
    )
  }

  sel$tiff_file[[1]]
}

.rects_from_df <- function(df) {
  lapply(seq_len(nrow(df)), function(i) {
    TiffPlotR::TiffRect(df$XMin[i], df$XMax[i], df$YMin[i], df$YMax[i])
  })
}

#' Load Selected Cell-Level Data For a CellQuery
#'
#' Reads cell-level CSV files referenced by selected samples in a
#' [CellQueryInfo-class] object and combines them into one data frame.
#'
#' @param object A [CellQueryInfo-class] object.
#'
#' @return Data frame containing bound rows from selected cell data files,
#'   including a `sample_id` column.
#' @examples
#' \dontrun{
#' q <- CellQuery()
#' cell_df <- load_query_cell_data(q)
#' }
#' @export
load_query_cell_data <- function(object) {
    stopifnot(methods::is(object, "CellQueryInfo"))

    cell_info_df <- get_query_cell_files_df(object)
    if (!nrow(cell_info_df)) {
        return(data.frame())
    }

    to_load <- cell_info_df$file
    names(to_load) <- cell_info_df$sample_id

    cell_df_l <- lapply(to_load, function(f) {
        f <- file.path(get_wrangled_cell_data_dir(), f)
        df <- readr::read_csv(f, show_col_types = FALSE)
        if ("AlgorithmName" %in% colnames(df)) {
            df$AlgorithmName <- NULL
        }
        if ("...1" %in% colnames(df)) {
            df$...1 <- NULL
        }
        df
    })

    dplyr::bind_rows(cell_df_l, .id = "sample_id")
}

#' Select Representative Cells Per Sample
#'
#' Filters cells by a marker classification/value and samples up to `n_cells`
#' rows per sample.
#'
#' @param cell_df Cell-level data frame or a [CellDataStore-class] object.
#' @param marker_col Column name used for filtering.
#' @param marker_value Value in `marker_col` to keep.
#' @param n_cells Maximum number of cells sampled per sample.
#' @param sample_col Column name for sample id grouping.
#' @param seed Integer random seed used for reproducible sampling.
#'
#' @return A named list of data frames split by sample id.
#' @examples
#' \dontrun{
#' q <- CellQuery()
#' cell_df <- load_query_cell_data(q)
#' reps <- select_representative_cells(cell_df)
#' }
#' @export
select_representative_cells <- function(
        cell_df,
        marker_col = "Opal520Classification",
        marker_value = 1,
        n_cells = 9,
        sample_col = "sample_id",
        seed = 1
) {
    if (methods::is(cell_df, "CellDataStore")) {
        cell_df <- get_selected_cell_data(cell_df)
    }

    stopifnot(marker_col %in% colnames(cell_df))
    stopifnot(sample_col %in% colnames(cell_df))

    query_df <- cell_df %>% dplyr::filter(.data[[marker_col]] == marker_value)
    query_df_l <- split(query_df, query_df[[sample_col]])

    set.seed(seed)
    query_df_l <- lapply(query_df_l, function(x) {
        if (!nrow(x)) {
            return(x)
        }
        n_take <- min(nrow(x), n_cells)
        x[sample(nrow(x), size = n_take), ]
    })

    query_df_l
}

#' Fetch Representative TIFF Crops For Selected Cells
#'
#' Uses selected sample TIFF paths in a [CellQueryInfo-class] object and
#' extracts representative annotated cell crops for each sample.
#'
#' @param object A [CellQueryInfo-class] object.
#' @param sampled_cells Either a data frame containing a `sample_id` column or
#'   a named list of data frames split by sample id.
#' @param fetch_resize_mult Expansion multiplier for each cell bounding box
#'   before TIFF extraction.
#' @param max_images_per_sample Maximum number of representative cells/TIFF crops
#'   to return per sample.
#' @param annotate_color Annotation color passed to
#'   [TiffPlotR::rect_annotate()].
#'
#' @return A named list keyed by sample id where each element is a list of
#'   `TiffPlotR` image objects.
#' @examples
#' \dontrun{
#' q <- CellQuery()
#' cell_df <- load_query_cell_data(q)
#' reps <- select_representative_cells(cell_df)
#' imgs <- fetch_representative_tiff_images(q, reps)
#' plot(imgs[[1]][[1]])
#' }
#' @export
fetch_representative_tiff_images <- function(
        object,
        sampled_cells,
        fetch_resize_mult = 2,
        max_images_per_sample = 3,
        red_channel = 2,
        green_channel = 3,
        blue_channel = 1,
        annotate_color = "yellow"
) {
    stopifnot(methods::is(object, "CellQueryInfo"))
    .stop_if_missing_tiffplotr()

    if (is.data.frame(sampled_cells)) {
        stopifnot("sample_id" %in% colnames(sampled_cells))
        sampled_cells <- split(sampled_cells, sampled_cells$sample_id)
    }

    sample_ids <- names(sampled_cells)
    image_files <- setNames(
        sapply(sample_ids, function(x) .find_tiff_file_by_sample(object, x)),
        sample_ids
    )

    image_res <- lapply(sample_ids, function(sample_id) {
        sample_df <- sampled_cells[[sample_id]]
        img_file <- image_files[[sample_id]]

        if (is.na(img_file) || !nrow(sample_df)) {
            return(list())
        }

        rects <- .rects_from_df(sample_df)
        rects <- rects[seq_len(min(length(rects), max_images_per_sample))]

        max_dim <- max(vapply(rects, function(r) {
            max(r@xmax - r@xmin, r@ymax - r@ymin)
        }, numeric(1)))
        resize_dim <- max_dim * fetch_resize_mult

        lapply(rects, function(r) {
            r_fetch <- TiffPlotR::rect_resize_abs(r, resize_dim, resize_dim)
            img_res <- TiffPlotR::fetchTiffData.rgb(
                img_file,
                r_fetch,
                blue_channel = blue_channel,
                red_channel = red_channel,
                green_channel = green_channel
            )
            img_res@plots$rgb <- TiffPlotR::rect_annotate(
                img_res@plots$rgb,
                r,
                color = annotate_color
            )
            img_res
        })
    })
    names(image_res) <- sample_ids

    image_res
}
