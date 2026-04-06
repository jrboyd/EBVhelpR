library(TiffPlotR)
library(EBVhelpR)
library(tidyverse)

#### setup ####
res_dir = "output_cell_views_04092026_v1"
res_file = function(f){
    f = file.path(dirname(f), paste0(plot_group, "_", f))
    out_f = file.path(res_dir, f)
    dir.create(dirname(out_f), recursive = TRUE, showWarnings = FALSE)
    out_f
}
#plots will be prefixed by ## id for organization
plot_group = "00"
increase_plot_group = function(){
    str = format(as.numeric(plot_group) + 1, digits = 2, width = 2)
    str = gsub(" ", "0", str)
    plot_group <<- str
}

check_saveplot = function(name){
    file.exists(res_file(paste0(name, ".png")))
}

my_saveplot = function(plot, name, width, height){
    ggsave(res_file(paste0(name, ".png")), plot, width = width, height = height, bg = "white")
    # ggsave(res_file(paste0(name, ".pdf")), plot, width = width, height = height)
}

#### load ####
cq_r4 = readRDS("output_presentation_04092026_v2/05_cq_r4.Rds")
cq_r3i = readRDS("output_presentation_04092026_v2/05_cq_r3i.Rds")
cq_r3i@summary_df %>% head
cq_r3i@all_cell_files_df$file

cq_p = readRDS("output_presentation_04092026_v2/05_cq_p.Rds")
get_tiff_file_path_df(cq_p)
tmp = get_query_tiff_paths_df(cq_p)
TiffPlotR::read_tiff_meta_data(tmp$tiff_file[1])
get_tiff_file_path_df(tmp$tiff_file[1])

cq = cq_r4
#update tiff paths to locally available
local_tiff_df = EBVhelpR:::get_tiff_file_path_df()
local_tiff_df = EBVhelpR:::.df_prep(local_tiff_df)
cq@tiff_paths_df = local_tiff_df
cq@tiff_paths_df = cq@tiff_paths_df %>% filter(unique_id %in% cq@selected_unique_ids) %>% filter(assay %in% cq@assay_type)
cq@tiff_paths_df = cq@tiff_paths_df[!duplicated(cq@tiff_paths_df$unique_id),]

cq@all_cell_files_df = cq@all_cell_files_df[!duplicated(cq@all_cell_files_df$unique_id),]
# cq@tiff_paths_df[duplicated(cq@tiff_paths_df$unique_id),]

assay_name = cq@assay_type



cq = filter_query_to_tiff_path_samples(cq)

qsum = get_query_summary_df(cq)
qcell_files = get_query_cell_files_df(cq)

todo = qcell_files$unique_id %>% unique

sel_id = todo[1]


for(sel_id in todo){
    message(sel_id)
    if(check_saveplot(paste0("validate_", sel_id))){
        message("skipping")
        next
    }

    cq.i = set_selected_unique_ids(cq, sel_id)
    cell_df = load_query_cell_data(cq.i)
    # qcell_files = get_query_cell_files_df(cq.i)
    tiff_df = get_query_tiff_paths_df(cq.i)
    tiff_f = tiff_df$tiff_file[1]
    stopifnot(length(tiff_f) == 1)

    cell_df = cell_df %>% mutate(x = (XMin + XMax)/2, y = (YMin + YMax)/2)


    img_info = read_tiff_meta_data(tiff_f)
    max_x = max(img_info$sizeX)
    max_y = max(img_info$sizeY)
    min_dim = min(max_x, max_y)

    # debug(fetchTiffData.rgb)
    img_res = tryCatch({
        fetchTiffData.rgb(tiff_path = tiff_f, channel_names = EBV_CHANNELS[[assay_name]],
                                    blue_channel = 1, red_channel = 6, green_channel = 2)

    }, error = function(e){
        NULL
    })
    if(is.null(img_res)){
        message("error during loading.")
        next
    }

    tp = sample(nrow(cell_df), 1e3)
    p1 = img_res@plots$rgb
    p2 = img_res@plots$rgb +
        annotate("point", x= cell_df$x[tp], y = cell_df$y[tp], color = "magenta") +
        coord_fixed()
    pg_validate = cowplot::plot_grid(p1, p2)

    eber_status = get_query_summary_df(cq.i)$EBER_status[1]
    p_title = ggplot() + labs(title = paste(sel_id, ":", eber_status)) +
        theme_void()
    pg_final = cowplot::plot_grid(p_title, pg_validate, ncol = 1, rel_heights = c(1, 12))
    my_saveplot(pg_final, paste0("validate_", sel_id), width = 9.8, 5.2)
    message("saved plot")
    # plot(pg_final)
}

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
#
# nbin = starting_nbin
# bin_size_x = diff(range(cell_df$x)) / nbin
# bin_size_y = diff(range(cell_df$y)) / nbin
# bin_size = min(bin_size_x, bin_size_y)
# cell_df = cell_df %>% mutate(xbin = round(x / bin_size), ybin = round(y / bin_size))
# bin_sizes = cell_df %>% group_by(xbin, ybin) %>% summarise(N = length(xbin)) %>% arrange(-N)
# bin_sizes
# bin_factor = sqrt(max(bin_sizes$N) / cell_target)
# nbin = nbin * bin_factor
#
# bin_size_x = diff(range(cell_df$x)) / nbin
# bin_size_y = diff(range(cell_df$y)) / nbin
# bin_size = min(bin_size_x, bin_size_y)
# cell_df = cell_df %>% mutate(xbin = round(x / bin_size), ybin = round(y / bin_size))
# bin_sizes = cell_df %>% group_by(xbin, ybin) %>% summarise(N = length(xbin)) %>% arrange(-N)
# bin_sizes
#
# sel_cell_df = cell_df %>% filter(xbin == bin_sizes$xbin[1] & ybin == bin_sizes$ybin[1])
cq.i = set_selected_unique_ids(cq, sel_id)
cell_df = load_query_cell_data(cq.i)
# qcell_files = get_query_cell_files_df(cq.i)
tiff_df = get_query_tiff_paths_df(cq.i)
tiff_f = tiff_df$tiff_file[1]
stopifnot(length(tiff_f) == 1)

cell_df = cell_df %>% mutate(x = (XMin + XMax)/2, y = (YMin + YMax)/2)


sel_cell_df = locate_cell_view_by_density(cell_df, n_cells = 10, bin_size = 100)

r = TiffRect(min(sel_cell_df$XMin), max(sel_cell_df$XMax), min(sel_cell_df$YMin), max(sel_cell_df$YMax))
w = r@coords$xmax - r@coords$xmin
h = r@coords$ymax - r@coords$ymin
#make rect square
r = rect_resize_abs(r, max(w, h), max(w, h))
r = rect_resize_abs(r, 500, 500)


fetchTiffData.rgb(tiff_path = tiff_f, rect = r,
                  channel_names = EBV_CHANNELS[[assay_name]],
                  blue_channel = 1, red_channel = 6, green_channel = 3,
                  resolution = 1)

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
img_res.z@plots$rgb %>% rect_annotate(view_rect, color = r_color) +
    coord_fixed(xlim = c(r@coords$xmin, r@coords$xmax),
                ylim = c(r@coords$ymin, r@coords$ymax))

#### rename channels
opal_codes = EBV_OPAL_DECODE[[cq.i@assay_type]]
for(name in names(opal_codes)){
    colnames(cell_df) = sub(name, opal_codes[name], colnames(cell_df))
}
colnames(cell_df)

cell_df.eber_pos = filter(cell_df, EBERClassification == 1)
sel_cell_df = locate_cell_view_by_density(cell_df.eber_pos, n_cells = 3, bin_size = 100)

#### load at high density area for signal ####

tp = sample(nrow(cell_df), 1e3)
img_res.z@plots$rgb +
    # annotate("point", x= cell_df$x[tp], y = cell_df$y[tp], color = "green") +
    coord_fixed(xlim = c(r@coords$xmin, r@coords$xmax), ylim = c(r@coords$ymin, r@coords$ymax))

r2 = rect_resize_mult(r, .2)
fetchTiffData.rgb(tiff_path = tiff_f, rect = r2,
                  channel_names = EBV_CHANNELS[[assay_name]], blue_channel = 1, red_channel = 4, green_channel = 5)

TiffRect()
