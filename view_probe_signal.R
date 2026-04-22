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
cq_p@tiff_paths_df

get_query_tiff_paths_df(cq_p)
get_query_cell_files_df(cq_p)
cq_4 = CellQuery(assay_type = EBV_ASSAY_TYPES$RNAScope_4plex) %>% filter_query_to_tiff_path_samples()
get_query_tiff_paths_df(cq_4)
cq_3i = CellQuery(assay_type = EBV_ASSAY_TYPES$`RNAScope_3plex+IF`) %>% filter_query_to_tiff_path_samples()
get_query_tiff_paths_df(cq_3i)

# limit to validated cell alignments
valid_cell_pngs = dir("validated_cells/", recursive = TRUE, pattern = ".png$", full.names = TRUE)
valid_df = data.frame(valid_file = valid_cell_pngs)
valid_df = valid_df %>% mutate(assay = basename(dirname(valid_file)) %>% sub("_JRB", "", .), unique_id = basename(valid_file) %>% sub("00_validate_", "", .) %>% sub(".png", "", .))
valid_ids.by_assay = split(valid_df$unique_id, valid_df$assay)
valid_ids.by_assay$Phenocycler

all_cqs = list(cq_4, cq_3i, cq_p)
names(all_cqs) = sapply(all_cqs, function(cq)cq@assay_type)

stopifnot(setequal(names(all_cqs), names(valid_ids.by_assay)))

cq = all_cqs$Phenocycler
all_cqs = lapply(all_cqs, function(cq){
    cq = set_selected_unique_ids(cq, valid_ids.by_assay[[cq@assay_type]])
    cq = filter_query_to_tiff_path_samples(cq)
    cq
})

lapply(all_cqs, get_query_cell_files_df)

cq = all_cqs$`RNAScope_3plex+IF`

bin_cell_df = function(cell_df, bin_size = 100){
    cell_df = cell_df %>% mutate(x = (XMin + XMax)/2, y = (YMin + YMax)/2)
    cell_df = cell_df %>% mutate(xbin = round(x / bin_size), ybin = round(y / bin_size))
    cell_df
}

get_density_df = function(cell_df, bin_size = 100){
    cell_df = bin_cell_df(cell_df, bin_size)
    bin_sizes = cell_df %>% group_by(xbin, ybin) %>% summarise(N = length(xbin)) %>% arrange(-N)
    bin_sizes
}

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
locate_cell_view_by_density = function(cell_df, n_cells = 10, bin_size = 100, n_views = 1){
    cell_df = bin_cell_df(cell_df, bin_size)
    bin_sizes = cell_df %>% group_by(xbin, ybin) %>% summarise(N = length(xbin)) %>% arrange(-N)
    rng = range(bin_sizes$N)
    if(n_cells < min(rng) | n_cells > max(rng)){
        warning("cell density not obtainable. observed range is ", min(rng), " to ", max(rng))
    }
    bin_sizes = bin_sizes %>% arrange(abs(N - n_cells))
    bin_sizes = bin_sizes %>% mutate(view_id = paste(xbin, ybin))
    bin_sizes.sel = bin_sizes[seq(1, n_views),]
    cell_df.sel = merge(cell_df, bin_sizes.sel)
    # sel_cell_df = cell_df %>% filter(xbin == bin_sizes$xbin[1] & ybin == bin_sizes$ybin[1])
    cell_df.sel
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
    if(!is.null(rownames(x))){
        rect_x$name = rownames(x)
    }else{
        is_not_rect = !colnames(x) %in% match_in_x
        name_col = colnames(x)[is_not_rect][1]
        rect_x$name = x[[name_col]]
    }
    do.call(TiffPlotR::TiffRect, rect_x)
}

rect_range = function(r){
    TiffRect(min(r@coords$xmin), max(r@coords$xmax),
             min(r@coords$ymin), max(r@coords$ymax))
}

library(TiffPlotR)

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
        img_data = TiffPlotR::fetchTiffData(tiff_f, rect = view_rects[[sel_id]], channel_names = chan_names, selected_channels = chan_names)



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

cq = all_cqs$RNAScope_4plex
all_vplots = lapply(all_cqs, master_plot_validation)
all_vplots$RNAScope_4plex$CTEBV_11

view_cell_density = function(cq, bin_sizes){
    cdata = load_query_cell_data(cq)
    # cdata %>% select(contains("Positive")) %>% select(contains("LMP"))
    # cdata = cdata %>% filter(LMP1.Positive.Classification == 1)
    cdata.by_id = split(cdata, cdata$unique_id)
    density_df.by_id = lapply(cdata.by_id, get_density_df)
    plot_df = bind_rows(density_df.by_id, .id = "sample_id")
    ggplot(plot_df, aes(x = N, fill = sample_id)) +
        geom_histogram(binwidth = 1) +
        facet_grid(sample_id~., scales = "free_y")
}

view_cell_density.positive = function(cq, bin_sizes){
    cdata = load_query_cell_data(cq)
    all_cn = colnames(cdata)
    all_cn = all_cn[grepl("Positiv", all_cn)]
    names(all_cn) = all_cn
    all_plots = lapply(all_cn, function(cn){
        cdata.sel = cdata %>% filter(!!sym(cn) == 1)
        cdata.by_id = split(cdata.sel, cdata.sel$unique_id)
        density_df.by_id = lapply(cdata.by_id, get_density_df)
        plot_df = bind_rows(density_df.by_id, .id = "sample_id")
        plot_df = merge(plot_df, cq@meta_data_df)
        ggplot(plot_df, aes(x = N, fill = EBER_status)) +
            geom_histogram(binwidth = 1) +
            facet_grid(sample_id~., scales = "free_y") +
            scale_fill_manual(values = EBVhelpR::get_colors_EBER_status()) +
            labs(title = cn)
    })
    all_plots
}

region_colors = seqsetvis::safeBrew(names(pos_cns.by_group))
region_colors["Classification"] = "yellow"
region_colors["Cytoplasm"] = "green"
region_colors["Nucleus"] = "cyan"


region_sizes = region_colors
region_sizes["Classification"] = 1
region_sizes["Cytoplasm"] = .9
region_sizes["Nucleus"] = 1/.9

master_plot_density = function(cq, target_density, n_views){
    cdata = load_query_cell_data(cq)
    cns = colnames(cdata)
    pos_cns = cns[grepl("Posit", cns)]
    groups = sapply(strsplit(sub("c.Myc", "CMyc", pos_cns), "\\."), function(x)x[3])
    pos_cns.by_group = split(pos_cns, groups)


    cdata.by_id = split(cdata, cdata$unique_id)
    cdata_sel.by_id = lapply(cdata.by_id, function(cell_data){
        locate_cell_view_by_density(cell_data, n_cells = 20, bin_size = 100)
    })
    # rects.by_id$D_EB_18$XMin
    # rects.by_id$D_EB_18$XMax
    # rects.by_id$D_EB_18$YMin
    # rects.by_id$D_EB_18$YMax
    # rects.by_id[[sel_id]]
    rects.by_id = lapply(cdata_sel.by_id, as.TiffRect)
    view_rects = lapply(rects.by_id, rect_range)
    view_rects$D_EB_18


    todo = names(view_rects)
    names(todo) = todo
    sel_id = todo[1]

    validation_img_dat = lapply(todo, function(sel_id){
        view_r = view_rects[[sel_id]]
        tiff_paths_df.sel = tiff_paths_df %>% dplyr::filter(unique_id == sel_id)
        stopifnot(nrow(tiff_paths_df.sel) == 1)
        tiff_f = tiff_paths_df.sel$tiff_file[1]
        chan_names = EBVhelpR::EBV_CHANNELS[[cq@assay_type]]

        img_data = TiffPlotR::fetchTiffData(tiff_f, rect = view_r, channel_names = chan_names, selected_channels = chan_names)


        # ggplot() %>% rect_annotate(view_r) %>% rect_annotate(rects.by_id[[sel_id]], color = "black")
        cell_data = cdata.by_id[[sel_id]] %>% as.data.frame

        rownames(cell_data) = paste0("cell_", seq_len(nrow(cell_data)))
        cell_data %>% as.TiffRect()
        ann_cells = cell_data %>% as.TiffRect() %>% TiffPlotR::rect_test_overlap(., view_r, subset = TRUE)
        cell_data.sel = cell_data[ann_cells@coords$name,]

        img_data@plots$normalized.annotated = rect_annotate(img_data@plots$normalized, ann_cells, color =  "lightgreen", alpha = .3)  + labs(title = sel_id)
        cell_data.sel$i = 1
        cell_data.sel$j = 1
        cell_data.sel$norm_value = 1
        cell_data.sel = cell_data.sel %>% select(XMin, XMax, YMin, YMax, all_of(pos_cns)) %>% pivot_longer(all_of(pos_cns))
        cell_data.sel = subset(cell_data.sel, value == 1)
        call_rects = cell_data.sel %>%
            mutate(name = sub("c.Myc", "c-Myc", name)) %>%
            mutate(name = sub("Pax5", "PAX5", name)) %>%
            mutate(channel = sub(".Pos.+", "", name)) %>%
            mutate(region = stringr::str_split_i(name, "\\.", 3))
        call_rects
        if(!all(call_rects$channel %in% chan_names)){
            setdiff(call_rects$channel, chan_names)
            browser()
        }
        call_rects$channel = factor(call_rects$channel, levels = chan_names)

        call_rects$region %>% table
        split(call_rects, call_rects$region)

        call_rects$i = 1
        call_rects$j = 1
        call_rects$norm_value = 1
        call_rects = as.data.frame(call_rects)
        rownames(call_rects) = paste0("rect_", seq_len(nrow(call_rects)))
        call_rects.by_region = split(call_rects, call_rects$region)
        x = call_rects.by_region[[1]]
        for(region_name in names(call_rects.by_region)){
            x = call_rects.by_region[[region_name]]
            new_r = rect_resize_mult(as.TiffRect(x), region_sizes[region_name])
            x$XMin = new_r@coords$xmin
            x$XMax = new_r@coords$xmax
            x$YMin = new_r@coords$ymin
            x$YMax = new_r@coords$ymax
            call_rects.by_region[[region_name]] = x
        }
        call_rects.plot = bind_rows(call_rects.by_region)
        img_data@plots$normalized.calls = img_data@plots$normalized +
            geom_rect(data = call_rects.plot, aes(xmin = XMin, xmax = XMax, ymin = YMin, ymax = YMax, color = region), fill = NA) +
            scale_color_manual(values = region_colors)
        img_data@activePlot = "normalized.calls"
        # plot(img_data)
        img_data
    })

    all_cell_plots = lapply(validation_img_dat, function(img_data){
        img_data@plots$normalized.annotated + labs(caption = "annotated with cell segmentation")
    })

    all_probe_plots = lapply(validation_img_dat, function(img_data){
        img_data@plots$normalized.calls + labs(caption = "annotated with probe classification")
    })
    all_probe_plots$D_EB_18
    all_probe_plots$D_EB_20
    all_probe_plots$D_EB_8
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
