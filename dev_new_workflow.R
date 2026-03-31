# library(magrittr)
# library(tidyverse)
library(EBVhelpR)

# Use live package sources in development so newly added functions/classes are available.
if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE)
}

if (!requireNamespace("TiffPlotR", quietly = TRUE)) {
  stop("Package `TiffPlotR` is required for TIFF image extraction.", call. = FALSE)
}

# Resolve package symbols even if they are not yet exported in the installed build.
.get_pkg_fn <- function(name) {
  ns <- asNamespace("EBVhelpR")
  get(name, envir = ns)
}

# CellQuery <- .get_pkg_fn("CellQuery")
# set_selected_sample_ids <- .get_pkg_fn("set_selected_sample_ids")
# get_query_summary_df <- .get_pkg_fn("get_query_summary_df")
# load_query_cell_data_store <- .get_pkg_fn("load_query_cell_data_store")
# get_full_cell_data <- .get_pkg_fn("get_full_cell_data")
# get_selected_cell_data <- .get_pkg_fn("get_selected_cell_data")
# get_original_cells_per_sample_id <- .get_pkg_fn("get_original_cells_per_sample_id")
# restore_selected_cell_data <- .get_pkg_fn("restore_selected_cell_data")
# replace_full_with_selected_cell_data <- .get_pkg_fn("replace_full_with_selected_cell_data")
# select_representative_cells <- .get_pkg_fn("select_representative_cells")
# fetch_representative_tiff_images <- .get_pkg_fn("fetch_representative_tiff_images")

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
  filter(Opal520Classification == 1) %>%
  mutate(any_marker = Opal520Classification + Opal620Classification)

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

debug(fetch_representative_tiff_images)
img_res.with_bg <- fetch_representative_tiff_images(
    object = cq,
    sampled_cells = reps,
    fetch_resize_mult = 2,
    max_images_per_sample = 3,
    annotate_color = "yellow",
    cell_data_store = cds
)

img_res.with_bg$CTEBV_11[[2]]

library(ggplot2)

plots.l = lapply(img_res, function(tiff_s){
    lapply(tiff_s, function(x){
        p = x@plots$rgb
        p + theme_void() + guides(color = "none")
    })
})
plots = unlist(plots.l)
cowplot::plot_grid(plotlist = plots, ncol = 3)

