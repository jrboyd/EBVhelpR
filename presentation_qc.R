library(EBVhelpR)
library(tidyverse)

# undebug(load_cell_source_files)
cell_sources = load_cell_source_files()
sel_ids = readRDS("output_presentation_04092026_v6/run_all_samples/03_final_ids.Rds")

cell_sources = cell_sources %>% filter(unique_id %in% sel_ids)
# ddir = get_wrangled_cell_data_dir()
#
# todo = cell_sources$file
# names(todo) = paste(cell_sources$assay, cell_sources$unique_id)
# todo = file.path(ddir, todo)
#
# count_file = dir(ddir, pattern = "cell_counts.txt", full.names = TRUE)
# tmp = read.table(count_file, sep = "\n")
# tmp$V1 = sub(" +", "", tmp$V1)
# tmp$V1 = sub(" ", "\t", tmp$V1)
# cell_counts = tmp %>% separate(V1, c("cell_count", "file"), sep = "\t")
# cell_counts$cell_count = as.numeric(cell_counts$cell_count)

# cell_counts = merge(cell_sources, cell_counts)
cell_sources
cell_o = cell_sources %>% group_by(unique_id) %>% summarize(cell_count = mean(cell_count)) %>%
    arrange(-cell_count)

cell_sources$unique_id = factor(cell_sources$unique_id, levels = cell_o$unique_id)

cell_sources

ggplot(cell_sources, aes(x = unique_id, y = cell_count, fill = sample_type)) +
    geom_col() +
    facet_wrap(~assay, ncol = 1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8))

tiff_paths = EBVhelpR::get_tiff_file_path_df()
comb_df = merge(cell_sources, tiff_paths)

low_cell_tiff_df = comb_df %>% filter(sample_id == "D_EB_2")



cq_p = CellQuery(assay_type = EBV_ASSAY_TYPES$Phenocycler) %>% filter_query_to_tiff_path_samples()
get_query_tiff_paths_df(cq_p)
get_query_cell_files_df(cq_p)
cq_4 = CellQuery(assay_type = EBV_ASSAY_TYPES$RNAScope_4plex) %>% filter_query_to_tiff_path_samples()
get_query_tiff_paths_df(cq_4)
cq_3i = CellQuery(assay_type = EBV_ASSAY_TYPES$`RNAScope_3plex+IF`) %>% filter_query_to_tiff_path_samples()
get_query_tiff_paths_df(cq_3i)

cq = cq_p

#' Title
#'
#' @param cell_df input cell_df
#' @param n_cells number of cells
#' @param bin_size region area is bin_sizeXbin_size pixels
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

as.TiffRect = function(x){
    rect_cn = c("xmin", "xmax", "ymin", "ymax")
    match_in_x = sapply(rect_cn, function(cn){
        match_i = which(cn == tolower(colnames(x)))
        colnames(x)[match_i]
    })
    tmp = names(match_in_x)
    names(tmp) = match_in_x

    rect_x = x[, match_in_x]
    colnames(rect_x) = tmp[colnames(rect_x)]
    if(is.null(rownames(x))){
        rect_x$name = rownames(x)
    }else{
        is_not_rect = !colnames(x) %in% match_in_x
        name_col = colnames(x)[is_not_rect][1]
        rect_x$name = x[[name_col]]
    }
    do.call(TiffPlotR::TiffRect, rect_x)
}

library(TiffPlotR)

cq = cq_p

#validates match between cell data and images
#overlays n_cells
# focuses on region in image where cells are located
master_plot_validation = function(cq){
    tiff_paths_df = get_query_tiff_paths_df(cq)
    cdata = load_query_cell_data(cq)
    view_rects = cdata %>% group_by(unique_id) %>% summarise(XMin = min(XMin), XMax = max(XMax), YMin = min(YMin), YMax = max(YMax))
    view_rects = lapply(split(view_rects, view_rects$unique_id), as.TiffRect)

    sel_id = names(view_rects)[1]
    todo = names(view_rects)
    names(todo) = todo
    validation_plots = lapply(todo, function(sel_id){
        tiff_paths_df.sel = tiff_paths_df %>% dplyr::filter(unique_id == sel_id)
        stopifnot(nrow(tiff_paths_df.sel) == 1)
        tiff_f = tiff_paths_df.sel$tiff_file[1]
        chan_names = EBVhelpR::EBV_CHANNELS[[cq@assay_type]]
        img_data = TiffPlotR::fetchTiffData(tiff_f, rect = view_rects[[sel_id]], channel_names = chan_names)



        n_cells = 5e3
        safe_sample_n = function(tbl, size, ...){
            if(size > nrow(tbl)){
                size = nrow(tbl)
            }
            sample_n(tbl, size, ...)
        }
        ann_cells = cdata %>% filter(unique_id == sel_id) %>% safe_sample_n(n_cells) %>% as.TiffRect()
        img_data.ann = rect_annotate(img_data, ann_cells, annotate_center = TRUE, size = .2, alpha = .1)
        img_data.ann@plots$normalized.annotated + labs(title = sel_id)
    })
}

debug(master_plot_validation)
vplots_3i = master_plot_validation(cq_3i)
vplots_4 = master_plot_validation(cq_4)
vplots_p = master_plot_validation(cq_p)

TiffPlotR::fetchTiffData(tiff_f, rect = TiffRect(1.3e4, 1.8e4, 1.3e4, 1.8e4), channel_names = chan_names)
TiffPlotR::fetchTiffData(tiff_f, rect = TiffRect(1.3e4, 1.4e4, 1.3e4, 1.4e4), channel_names = chan_names, resolution = 1)

master_plot_density = function(cq, target_density, n_views){

}

# i = 2
# img_f = low_cell_tiff_df$tiff_file[i]
# cell_f = low_cell_tiff_df$file[i]
# cell_f = file.path(get_wrangled_cell_data_dir(), cell_f)

CellQuery()
load_query_cell_data()

cell_data = readr::read_csv(cell_f)
cell_data

basename(img_f)
TiffPlotR::fetchTiffData.rgb(img_f)
