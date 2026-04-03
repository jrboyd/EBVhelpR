library(magrittr)
library(tidyverse)
library(EBVhelpR)
#TODO resolve inconsistent behaviors
# sample_id everywhere
# source vs assay confusion
#### look at example cells ####
# debug(CellQuery)


cq = CellQuery(assay_type = EBV_ASSAY_TYPES$`RNAScope_3plex+IF`)
cq = CellQuery(assay_type = EBV_ASSAY_TYPES$RNAScope_4plex)
# debug(CellQuery)
cq = CellQuery(assay_type = EBV_ASSAY_TYPES$Phenocycler)


cq = CellQuery()
cq@all_cell_files_df = subset(cq@all_cell_files_df, file != "RNAScope_4plex/D-EB-3.cell_data.csv")
cq@summary_df$sample_id %>% table
cq@summary_df$Sample %>% table
tmp = cq@summary_df %>% subset(sample_id %in% cq@tiff_paths_df$sample_id)
tmp %>% select(unique_id, sample_id, probe_control) %>% unique

tmp$unique_id %>% table
tmp$unique_id %>% table


cq@tiff_paths_df
# debug(filter_query_to_tiff_path_samples)
cq = filter_query_to_tiff_path_samples(cq)
cq

set.seed(0)
sel_ids = cq@tiff_paths_df$unique_id %>% unique

# debug(set_selected_sample_ids)
cq = set_selected_unique_ids(cq, sel_ids)

# debug(get_query_cell_files_df)
cell_df = get_query_cell_files_df(cq)
cell_df$sample_id %>% table
cell_df

cq@selected_unique_ids
cq@all_cell_files_df$unique_id %>% table

# undebug(load_query_cell_data)
cell_df = load_query_cell_data(cq)
cell_df$unique_id %>% table
# sampled_cells = select_representative_cells(cell_df)

# undebug(select_representative_cells)
# Example usage
query_df.l <- select_representative_cells(
    cell_df,
    marker_col = "Opal520Classification",
    marker_value = 1,
    n_cells = 3,
    seed = 123
)

lapply(query_df.l, nrow)

# debug(fetch_representative_tiff_images)
# debug(fetch_representative_tiff_images)

lapply(cq@tiff_paths_df$tiff_file, TiffPlotR::read_tiff_meta_data)
# debug(TiffPlotR::fetchTiffData)
# lapply(cq@tiff_paths_df$tiff_file, TiffPlotR::fetchTiffData)

names(query_df.l)
# debug(fetch_representative_tiff_images)
tiff_sample <- fetch_representative_tiff_images(
    object = cq,
    sampled_cells = query_df.l,
    fetch_resize_mult = 1.5,
    max_images_per_sample = 3
)

cell_df$unique_id %>% table
cq@tiff_paths_df



names(tiff_sample)
tiff_sample[[1]][[2]]
tiff_sample[[1]][[1]]
tiff_sample[[1]][[3]]
tiff_sample[[2]][[1]]
tiff_sample[[2]][[2]]
tiff_sample[[3]][[2]]
tiff_sample[[4]][[2]]
names(tiff_sample)

library(ggplot2)
tiff_s = tiff_sample[[1]]
x = tiff_s[[1]]

p_ex = tiff_sample[[1]][[1]]@plots[[1]]

p_leg = cowplot::get_legend(p_ex + theme(legend.title = element_blank()))
class(p_leg)

plots.l = lapply(tiff_sample, function(tiff_s){
    lapply(tiff_s, function(x){
        p = x@plots$rgb.annotated
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
pg_body = cowplot::plot_grid(plotlist = plots, ncol = 3)

plots.row_names = lapply(names(plots.l), function(name){
    ggplot() + annotate("text", x= 0, y = 0, label = name) + theme_void()
})
pg_row_names = cowplot::plot_grid(plotlist = plots.row_names, ncol = 1)

pg_main = cowplot::plot_grid(pg_row_names, pg_body, nrow = 1, rel_widths = c(1, 5))
pg_main
# cowplot::plot_grid(pg_main, p_leg, ncol = 1)

# Example plot retrieval
#
# p = img_res@plots$normalized
# rect = r
# rect_df <- data.frame(xmin = rect@xmin, xmax = rect@xmax, ymin = rect@ymin, ymax = rect@ymax)
# p + ggplot2::geom_rect(data = rect_df, mapping = ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE, color = "red")
