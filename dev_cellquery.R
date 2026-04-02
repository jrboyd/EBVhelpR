library(magrittr)
library(tidyverse)
library(EBVhelpR)
#TODO resolve inconsistent behaviors
# sample_id everywhere
# source vs assay confusion
write_package_data_for_file
#### look at example cells ####
# debug(CellQuery)
cq = CellQuery()
cq = CellQuery(assay_type = EBV_ASSAY_TYPES$`rnascope_3plex+IF`)
cq = CellQuery(assay_type = EBV_ASSAY_TYPES$phenocycler)
cq = CellQuery(assay_type = EBV_ASSAY_TYPES$rnascope_4plex)

set.seed(0)
sel_ids = cq@tiff_paths_df$sample_id %>% unique %>% sample(4)

# debug(set_selected_sample_ids)
cq = set_selected_sample_ids(cq, sel_ids)

cell_df = get_query_cell_files_df(cq)
cell_df$sample_id %>% table
cell_df %>% head

cell_df = load_query_cell_data(cq)
cell_df$sample_id %>% table
# sampled_cells = select_representative_cells(cell_df)

# Example usage
query_df.l <- select_representative_cells(
    cell_df,
    marker_col = "Opal520Classification",
    marker_value = 1,
    n_cells = 9,
    seed = 123
)

lapply(query_df.l, nrow)

# debug(fetch_representative_tiff_images)
tiff_sample <- fetch_representative_tiff_images(
    object = cq,
    sampled_cells = query_df.l,
    fetch_resize_mult = 1.5,
    max_images_per_sample = 3
)

tiff_sample$D_EB_18[[2]]

x = tiff_sample$D_EB_18[[1]]

plots.l = lapply(tiff_sample, function(tiff_s){
    lapply(tiff_s, function(x){
        p = x@plots$rgb
        p + theme_void() + guides(color = "none")
    })
})
plots = unlist(plots.l)

# plots = lapply(tiff_sample$D_EB_18, function(x){
#     p = x@plots$rgb
#     p + theme_void()
# })

#ensure consistent color scale
#ensure consistent pixel scale
cowplot::plot_grid(plotlist = plots, ncol = 3)

lapply(plots, )

tiff_sample$D_EB_20

# Example plot retrieval
head(query_df.l$CTEBV_11)
tiff_sample$CTEBV_11[[1]]
p = plot(tiff_sample$CTEBV_11[[1]])
class(p)
p + coord_fixed()
#
# p = img_res@plots$normalized
# rect = r
# rect_df <- data.frame(xmin = rect@xmin, xmax = rect@xmax, ymin = rect@ymin, ymax = rect@ymax)
# p + ggplot2::geom_rect(data = rect_df, mapping = ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE, color = "red")
