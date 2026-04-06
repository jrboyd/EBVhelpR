library(TiffPlotR)
library(EBVhelpR)
library(tidyverse)

cq_r4 = readRDS("output_presentation_04092026_v2/05_cq_r4.Rds")
cq_r3i = readRDS("output_presentation_04092026_v2/05_cq_r3i.Rds")
cq_r3i@summary_df %>% head
cq_r3i@all_cell_files_df$file



cq = cq_r3i
#update tiff paths to locally available
local_tiff_df = EBVhelpR:::get_tiff_file_path_df()
local_tiff_df = EBVhelpR:::.df_prep(local_tiff_df)
cq@tiff_paths_df = local_tiff_df
cq@tiff_paths_df = cq@tiff_paths_df %>% filter(unique_id %in% cq@selected_unique_ids) %>% filter(assay %in% cq@assay_type)

assay_name = cq@assay_type



cq = filter_query_to_tiff_path_samples(cq)

qsum = get_query_summary_df(cq)
qcell_files = get_query_cell_files_df(cq)
qtiff = get_query_tiff_paths_df(cq)

todo = qcell_files$unique_id %>% unique
sel_id = todo[1]

cq.i = set_selected_unique_ids(cq, sel_id)
cell_df = load_query_cell_data(cq.i)
# qcell_files = get_query_cell_files_df(cq.i)
tiff_df = get_query_tiff_paths_df(cq.i)
tiff_f = tiff_df$tiff_file
stopifnot(length(tiff_f) == 1)

cell_df = cell_df %>% mutate(x = (XMin + XMax)/2, y = (YMin + YMax)/2)


img_info = read_tiff_meta_data(tiff_f)
max_x = max(img_info$sizeX)
max_y = max(img_info$sizeY)
min_dim = min(max_x, max_y)

# debug(fetchTiffData.rgb)
img_res = fetchTiffData.rgb(tiff_path = tiff_f, channel_names = EBV_CHANNELS[[assay_name]], blue_channel = 1, red_channel = 6, green_channel = 2)
tp = sample(nrow(cell_df), 1e3)
p1 = img_res@plots$rgb
p2 = img_res@plots$rgb +
    annotate("point", x= cell_df$x[tp], y = cell_df$y[tp], color = "magenta") +
    coord_fixed()
pg_validate = cowplot::plot_grid(p1, p2)
pg_validate

zr = TiffRect(0e4, 6e4, 2.5e3, 1e4)
img_res = fetchTiffData.rgb(tiff_path = tiff_f, channel_names = EBV_CHANNELS[[assay_name]], blue_channel = 1, red_channel = 6, green_channel = 2, rect = zr)
tp = sample(nrow(cell_df), 10e3)
p1 = img_res@plots$rgb +
    coord_fixed(xlim = c(zr@coords$xmin, zr@coords$xmax), ylim = c(zr@coords$ymin, zr@coords$ymax))
p2 = img_res@plots$rgb +
    annotate("point", x= cell_df$x[tp], y = cell_df$y[tp], color = "magenta") +
    coord_fixed(xlim = c(zr@coords$xmin, zr@coords$xmax), ylim = c(zr@coords$ymin, zr@coords$ymax))
pg_validate = cowplot::plot_grid(p1, p2)
pg_validate

#### look at high density area ####
cell_target = 50
starting_nbin = 50

n_cells = 10
bin_size = 100
#' Title
#'
#' @param cell_df input cell_df
#' @param n_cells number of cells
#' @param bin_size region area is bin_sizeXbin_size
#'
#' @returns
#' @export
#'
#' @examples
locate_cell_view_by_density = function(cell_df, n_cells = 10, bin_size = 100){
    cell_df = cell_df %>% mutate(xbin = round(x / bin_size), ybin = round(y / bin_size))
    bin_sizes = cell_df %>% group_by(xbin, ybin) %>% summarise(N = length(xbin)) %>% arrange(-N)
    rng = range(bin_sizes$N)
    if(n_cells < min(rng) | n_cells > max(rng)){
        stop("cell density not obtainable. observed range is ", min(rng), " to ", max(rng))
    }
    bin_sizes = bin_sizes %>% arrange(abs(N - n_cells))
    sel_cell_df = cell_df %>% filter(xbin == bin_sizes$xbin[1] & ybin == bin_sizes$ybin[1])
    sel_cell_df
}

nbin = starting_nbin
bin_size_x = diff(range(cell_df$x)) / nbin
bin_size_y = diff(range(cell_df$y)) / nbin
bin_size = min(bin_size_x, bin_size_y)
cell_df = cell_df %>% mutate(xbin = round(x / bin_size), ybin = round(y / bin_size))
bin_sizes = cell_df %>% group_by(xbin, ybin) %>% summarise(N = length(xbin)) %>% arrange(-N)
bin_sizes
bin_factor = sqrt(max(bin_sizes$N) / cell_target)
nbin = nbin * bin_factor

bin_size_x = diff(range(cell_df$x)) / nbin
bin_size_y = diff(range(cell_df$y)) / nbin
bin_size = min(bin_size_x, bin_size_y)
cell_df = cell_df %>% mutate(xbin = round(x / bin_size), ybin = round(y / bin_size))
bin_sizes = cell_df %>% group_by(xbin, ybin) %>% summarise(N = length(xbin)) %>% arrange(-N)
bin_sizes

sel_cell_df = cell_df %>% filter(xbin == bin_sizes$xbin[1] & ybin == bin_sizes$ybin[1])

sel_cell_df = locate_cell_view_by_density(cell_df, n_cells = 10, bin_size = 100)

r = TiffRect(min(sel_cell_df$XMin), max(sel_cell_df$XMax), min(sel_cell_df$YMin), max(sel_cell_df$YMax))
w = r@coords$xmax - r@coords$xmin
h = r@coords$ymax - r@coords$ymin
#make rect square
r = rect_resize_abs(r, max(w, h), max(w, h))
r = rect_resize_abs(r, 100, 100)

img_res.z = fetchTiffData.rgb(tiff_path = tiff_f, rect = r,
                              channel_names = EBV_CHANNELS[[assay_name]],
                              blue_channel = 1, red_channel = 2, green_channel = 3,
                              resolution = 1)

img_res.z
all_rect = TiffRect(cell_df$XMin, cell_df$XMax, cell_df$YMin, cell_df$YMax, cell_df$ObjectId)
view_rect = rect_test_overlap(all_rect, r, subset = TRUE)

r_color = seqsetvis::col2hex("green")
# r_color = paste0(r_color, "55")
random_colors = sample(colors(), view_rect@coords %>% nrow)
img_res.z@plots$rgb %>% rect_annotate(view_rect, color = random_colors) +
    coord_fixed(xlim = c(r@coords$xmin, r@coords$xmax),
                ylim = c(r@coords$ymin, r@coords$ymax))


#### load at high density area for signal ####

tp = sample(nrow(cell_df), 1e3)
img_res.z@plots$rgb +
    # annotate("point", x= cell_df$x[tp], y = cell_df$y[tp], color = "green") +
    coord_fixed(xlim = c(r@coords$xmin, r@coords$xmax), ylim = c(r@coords$ymin, r@coords$ymax))

r2 = rect_resize_mult(r, .2)
fetchTiffData.rgb(tiff_path = tiff_f, rect = r2,
                  channel_names = EBV_CHANNELS[[assay_name]], blue_channel = 1, red_channel = 4, green_channel = 5)

TiffRect()
