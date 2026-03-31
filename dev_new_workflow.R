# library(magrittr)
# library(tidyverse)
library(EBVhelpR)

if (!requireNamespace("TiffPlotR", quietly = TRUE)) {
  stop("Package `TiffPlotR` is required for TIFF image extraction.", call. = FALSE)
}

# 1) Build a CellQuery object for an assay.
cq <- CellQuery(assay_type = EBV_ASSAY_TYPES$rnascope_4plex)
print(cq)

# 2) Select a subset of samples to work on in this dev run.
set.seed(0)
sel_ids <- sample(unique(get_query_summary_df(cq)$sample_id), 4)
cq <- set_selected_sample_ids(cq, sel_ids)

# 3) Load selected cell files into a CellDataStore (full + selected views).
cds <- load_query_cell_data_store(cq)
print(cds)

# Original per-sample counts at load time.
orig_counts <- get_original_cells_per_sample_id(cds)
print(orig_counts)

# 4) dplyr verbs operate on selected view only.
cds <- cds %>%
  dplyr::filter(Opal520Classification == 1) %>%
    dplyr::mutate(any_marker = Opal520Classification + Opal620Classification)

cat("Rows after selected filtering:", nrow(get_selected_cell_data(cds)), "\n")
cat("Rows in full data remain:", nrow(get_full_cell_data(cds)), "\n")

# 5) Restore selected view back to full data.
cds <- restore_selected_cell_data(cds)
cat("Rows after restore:", nrow(get_selected_cell_data(cds)), "\n")

# 6) Re-filter and commit selected view into full data.
cds <- cds %>% dplyr::filter(Opal520Classification == 1)
cds <- replace_full_with_selected_cell_data(cds)
cat("Rows after commit full <- selected:", nrow(get_full_cell_data(cds)), "\n")

# 7) Sample representative cells and fetch representative TIFF crops.
reps <- select_representative_cells(
  cds,
  marker_col = "Opal520Classification",
  marker_value = 1,
  n_cells = 9,
  seed = 123
)

img_res <- fetch_representative_tiff_images(
  object = cq,
  sampled_cells = reps,
  fetch_resize_mult = 2,
  max_images_per_sample = 3,
  annotate_color = "yellow"
)

img_res$CTEBV_11[[2]]

# debug(fetch_representative_tiff_images)
img_res.with_bg <- fetch_representative_tiff_images(
    object = cq,
    sampled_cells = reps,
    fetch_resize_mult = 5,
    max_images_per_sample = 3,
    annotate_color = "yellow",
    cell_data_store = cds,
    background_cell_color = "orange"
)

cds.full <- load_query_cell_data_store(cq)
img_res.reanno = reannotate_image_list(img_res, cds.full, target_image_name = "test")
cds.full@selected_df %>% head
cds.full %>% dplyr::filter(Opal520Classification == 1)
cds.full %>% dplyr::filter(Opal570Classification == 1)
cds.full %>% dplyr::filter(Opal620Classification == 1)

img_res.reanno = reannotate_image_list(img_res,
                                       cds.full,
                                       annotation_color = "white",
                                       target_image_name = "test")
img_res.reanno = reannotate_image_list(img_res.reanno,
                                       cds.full %>% dplyr::filter(Opal520Classification == 1),
                                       annotation_color = "red",
                                       source_image_name = "test", target_image_name = "test")
img_res.reanno = reannotate_image_list(img_res.reanno,
                                       cds.full %>% dplyr::filter(Opal570Classification == 1),
                                       annotation_color = "green",
                                       source_image_name = "test", target_image_name = "test")

reannotate_image_list = function(img_res, cds, annotation_color = "white", source_image_name = "rgb", target_image_name = paste0(source_image_name, ".reann")){
  stopifnot(methods::is(cds, "CellDataStore"))
  stopifnot(requireNamespace("TiffPlotR", quietly = TRUE))

  full_df <- get_full_cell_data(cds)
  # Cache all cell rects per sample_id so repeated images from the same sample
  # do not rebuild the same TiffRect object over and over.
  rect_cache <- list()

  rects_from_df <- function(df) {
    TiffPlotR::TiffRect(
      df$XMin,
      df$XMax,
      df$YMin,
      df$YMax,
      name = as.character(df$ObjectId)
    )
  }

  get_fetch_rect <- function(img) {
    slot_names <- methods::slotNames(img)
    rect_slots <- c("fetch_rect", "fetchRect", "rect", "region", "crop_rect", "subset_rect", "tiff_rect")

    for (slot_name in rect_slots) {
      if (slot_name %in% slot_names) {
        slot_val <- methods::slot(img, slot_name)
        if (methods::is(slot_val, "TiffRect")) {
          return(slot_val)
        }
      }
    }

    # Fallback: derive the viewing window from plot limits.
    gb <- ggplot2::ggplot_build(img@plots[[source_image_name]])
    panel <- gb$layout$panel_params[[1]]

    x_range <- panel$x.range
    if (is.null(x_range)) {
      x_range <- panel$x_range
    }
    y_range <- panel$y.range
    if (is.null(y_range)) {
      y_range <- panel$y_range
    }

    if (is.null(x_range) || is.null(y_range)) {
      stop("Could not infer fetch rectangle from image object/plot.")
    }

    TiffPlotR::TiffRect(
      min(x_range),
      max(x_range),
      min(y_range),
      max(y_range),
      name = "focus"
    )
  }

  get_sample_rects <- function(sample_id) {
    key <- as.character(sample_id)
    if (!key %in% names(rect_cache)) {
      sample_df <- full_df[full_df$sample_id == key, , drop = FALSE]
      rect_cache[[key]] <<- rects_from_df(sample_df)
    }
    rect_cache[[key]]
  }

  annotate_one <- function(img, sample_id) {
    if (is.null(sample_id) || !nzchar(sample_id)) {
      stop("Could not infer sample_id for an image. Pass the full named image list from fetch_representative_tiff_images().")
    }

    fetch_rect <- get_fetch_rect(img)
    background_rects <- get_sample_rects(sample_id)
    # Keep only rectangles that intersect the current image window.
    background_rects <- TiffPlotR::rect_test_overlap(background_rects, fetch_rect, subset = TRUE)

    if (nrow(background_rects@coords) > 0) {
      primary_name <- fetch_rect@coords$name
      # Exclude the focal/primary rectangle name so only background cells remain.
      if (!is.null(primary_name) && length(primary_name) > 0) {
        coords_df <- background_rects@coords
        coords_df <- coords_df[
          !(coords_df[["name"]] %in% primary_name),
          ,
          drop = FALSE
        ]
        background_rects@coords <- coords_df
      }
    }

    p <- img@plots[[source_image_name]]
    if (nrow(background_rects@coords) > 0) {
      img@plots[[target_image_name]] <- TiffPlotR::rect_annotate(
        p,
        background_rects,
        color = annotation_color
      )
    } else {
      img@plots[[target_image_name]] <- p + ggplot2::labs(caption = "no annotations in view")
    }

    img@activePlot <- target_image_name
    img
  }

  walk_images <- function(x, sample_id = NULL) {
    if (methods::is(x, "TiffPlotData") && "plots" %in% methods::slotNames(x)) {
      return(annotate_one(x, sample_id = sample_id))
    }

    if (!is.list(x)) {
      return(x)
    }

    nms <- names(x)
    out <- lapply(seq_along(x), function(i) {
      x_i <- x[[i]]
      name_i <- if (!is.null(nms)) nms[[i]] else NULL

      next_sample_id <- sample_id
      # If sample_id was not passed down yet, infer it from named top-level nodes
      # (the list returned by fetch_representative_tiff_images is keyed by sample).
      if (is.null(next_sample_id) && !is.null(name_i) && nzchar(name_i) && name_i %in% full_df$sample_id) {
        next_sample_id <- name_i
      }

      walk_images(x_i, sample_id = next_sample_id)
    })
    if (!is.null(nms)) {
      names(out) <- nms
    }
    out
  }

  out <- walk_images(img_res)

  out
}

img_res.with_bg$CTEBV_11[[2]]

img_res.with_bg$CTEBV_11[[1]]


library(ggplot2)

plots.l = lapply(img_res.with_bg, function(tiff_s){
    lapply(tiff_s, function(x){
        p = x@plots$rgb
        p + theme_void() + guides(color = "none")
    })
})
plots = unlist(plots.l)
cowplot::plot_grid(plotlist = plots, ncol = 3)

