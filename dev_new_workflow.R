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

reannotate_image_list = function(img_res, cds, annotation_color = "yellow", source_image_name = "rgb", target_image_name = paste0(source_image_name, ".reann")){
    #iterate all img_res in list (or single img_res, then be sure to return single img_res, not list)
    #match sample in cds and find overlapping rects, annotate with annotation_color
    p = img_res$CTEBV_11[[1]]@plots[[source_image_name]]
    img_res$CTEBV_11[[1]]@plots[[target_image_name]] = p + rect_annotate()
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

