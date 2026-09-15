library(EBVhelpR)
library(TiffPlotR)
library(tidyverse)
# install.packages("tidyverse")
# ??replace_values
as.TiffRect = function(x, name_col = NULL){
    rect_cn = c("xmin", "xmax", "ymin", "ymax")
    match_in_x = sapply(rect_cn, function(cn){
        match_i = which(cn == tolower(colnames(x)))
        colnames(x)[match_i]
    })
    tmp = names(match_in_x)
    names(tmp) = match_in_x

    rect_x = x[, match_in_x]
    colnames(rect_x) = tmp[colnames(rect_x)]
    if(!is.null(rownames(x)) & !is(x, "tbl")){
        rect_x$name = rownames(x)
    }else{
        if(is.null(name_col)){
            is_not_rect = !colnames(x) %in% match_in_x
            name_col = colnames(x)[is_not_rect][1]
        }
        rect_x$name = x[[name_col]]
    }
    do.call(TiffPlotR::TiffRect, rect_x)
}


# cq_to_count = EBVhelpR::CellQuery(EBV_ASSAY_TYPES$RNAScope_4plex)


.count_in_best = function(
        cq,
        sel_id,
        cache_dir,
        block_size = 3e3,
        target_cell_count = 1e3
){
    cache_file = file.path(cache_dir, paste0(sel_id, ".", block_size, ".", target_cell_count, ".Rds"))
    if(file.exists(cache_file)){
        return(readRDS(cache_file))
    }else{

        cq_sel = EBVhelpR::set_selected_unique_ids(cq, sel_id)
        cds = load_query_cell_data(cq_sel)
        message(sel_id, " found ", nrow(cds), " cells")


        tiff_df = EBVhelpR::get_query_tiff_paths_df(cq_sel)

        decode_opal_colnames = function(cq_i, cell_df = NULL){
            if(is.null(cell_df)){
                cell_df = load_query_cell_data(cq_i)
            }
            opal_codes = EBV_OPAL_DECODE[[cq_i@assay_type]]
            for(name in names(opal_codes)){
                colnames(cell_df) = sub(name, paste0(opal_codes[name], "_"), colnames(cell_df))
            }
            colnames(cell_df)
            cell_df
        }

        cds = decode_opal_colnames(cq_sel, cds)
        if("Object.Id" %in% colnames(cds)){
            cds = cds %>% rename(ObjectId = Object.Id)
        }
        cds = cds %>% select(unique_id, ObjectId, XMin, XMax, YMin, YMax, contains("Classification"))
        #select only EBER
        # cds = cds %>% select(unique_id, ObjectId, XMin, XMax, YMin, YMax, contains("Classification"))


        quick_center = function(obj_data){
            obj_data %>% mutate(x = (XMin + XMax)/2, y = (YMin + YMax)/2, name = ObjectId, .keep = "none")
        }

        message("classifying cells...")
        centers_all = cds %>% quick_center()
        # debug(as.TiffRect)
        # centers_any = cds %>% filter(if_any(contains("Classification"), ~ (. == 1))) %>% quick_center
        #select only EBER
        centers_any = cds %>% filter(if_any(contains("EBER"), ~ (. == 1))) %>% quick_center
        class_cns = cds %>% select(contains("Classification")) %>% colnames
        names(class_cns) = sub("_Class.+", "", class_cns)
        cn = class_cns[1]
        centers_each = lapply(class_cns, function(cn){
            hit = cds %>% filter(!!sym(cn) == 1)
            if(nrow(hit) == 0) return(NULL)
            # hit %>% as.TiffRect(name_col = "ObjectId") %>% rect_centers() %>% tibble()
            quick_center(hit)
        })
        centers_todo = c(centers_each, list(all = centers_all, any = centers_any))
        drop = sapply(centers_todo, is.null)
        centers_todo = centers_todo[!drop]



        x_rng = centers_any$x %>% range
        y_rng = centers_any$y %>% range

        #view rects are half steps of black size
        x_breaks = seq(from = min(x_rng)-block_size/4, to = max(x_rng) + block_size/2, by = block_size/2)
        y_breaks = seq(from = min(y_rng)-block_size/4, to = max(y_rng) + block_size/2, by = block_size/2)
        x_lims = tibble(i = seq(length(x_breaks) -1)) %>% mutate(xmin = x_breaks[i], xmax = x_breaks[i+1]) %>% rename(xi = i)
        y_lims = tibble(i = seq(length(y_breaks) -1)) %>% mutate(ymin = y_breaks[i], ymax = y_breaks[i+1]) %>% rename(yi = i)

        count_rects = merge(x_lims, y_lims)
        message("counting cells by group...")
        for(name in names(centers_todo)){
            centers_i = centers_todo[[name]]
            lookup = "count" %>% setNames(paste0("count_", name))
            count_rects = count_rects %>% group_by(xi, yi) %>% mutate(count = sum(xmin < centers_i$x & xmax >= centers_i$x & ymin < centers_i$y & ymax >= centers_i$y)) %>% rename(all_of(lookup))
        }

        # count_rects


        # plot_df = pivot_longer(count_rects, starts_with("count"))
        #
        # ggplot(plot_df, aes(x = value)) +
        #     geom_histogram(bins = 200) +
        #     facet_grid(name~., scales = "free_y")

        #select count_rect with highest density
        count_rects = count_rects %>% mutate(density = count_any^2 / count_all)

        best_rect = count_rects %>% ungroup() %>% filter(density == max(density, na.rm = TRUE))

        #identify mean center of positive
        centers_in_best = lapply(centers_todo, function(c_df){
            c_df %>% filter(x > best_rect$xmin & x <= best_rect$xmax & y > best_rect$ymin & y <= best_rect$ymax)
        })

        center_point = centers_in_best$any %>% summarise(x = median(x), y = median(y))

        message("selecting nearest cells to center...")
        #select nearest target cells

        dist_df = centers_todo$all %>% mutate(distance = sqrt((center_point$x - x)^2 + (center_point$y - y)^2))
        dist_df = dist_df %>% mutate(rnk = rank(distance, ties.method = "first")) %>% arrange(rnk)
        best_cells_df = dist_df %>% filter(rnk <= target_cell_count)
        best_cell_ids = best_cells_df$name

        centers_nearest = lapply(centers_todo, function(x){
            x %>% filter(name %in% best_cell_ids)
        })

        # ch_res = chull(centers_nearest$all$x, centers_nearest$all$y)
        # #area of chull?
        # ch_poly = tibble(x = centers_nearest$all$x[ch_res],
        #                  y = centers_nearest$all$y[ch_res]
        # )

        hits_df = bind_rows(centers_nearest, .id = "type")
        count_df = hits_df %>% group_by(type) %>% summarise(N = length(type))
        count_df$unique_id = sel_id
        out = list(center = center_point, hits = hits_df, counts = count_df)
        saveRDS(out, cache_file)
        return(out)
    }
}



count_best_view_in_cq = function(cq_to_count){
    cache_dir = file.path(".cache_best_view.only_EBER", cq_to_count@assay_type)
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

    todo_ids = cq_to_count@selected_unique_ids
    names(todo_ids) = todo_ids

    count_res = lapply(todo_ids, .count_in_best, cq = cq_to_count, cache_dir = cache_dir)
    count_res
}


res_4plex = count_best_view_in_cq(EBVhelpR::CellQuery(EBV_ASSAY_TYPES$RNAScope_4plex))


# res_3plexIF = count_best_view_in_cq(EBVhelpR::CellQuery(EBV_ASSAY_TYPES$`RNAScope_3plex+IF`))
# res_P = count_best_view_in_cq(EBVhelpR::CellQuery(EBV_ASSAY_TYPES$Phenocycler))

id_file = "~/../OneDrive - UVM Larner College of Medicine/projects_ashley/EBV_DLBCL/paper1_enhanced_EBV_detection/output_presentation_04092026_v7/run_all_samples/03_final_ids.Rds"
if(!file.exists(id_file)){
    id_file = "~/../../OneDrive - UVM Larner College of Medicine/projects_ashley/EBV_DLBCL/paper1_enhanced_EBV_detection/output_presentation_04092026_v7/run_all_samples/03_final_ids.Rds"
}
final_ids = readRDS(id_file)

# final_ids = readRDS("../../projects_ashley/EBV_DLBCL/paper1_enhanced_EBV_detection/despiked_id_order.Rds")


wgs_files_df = EBVhelpR::setup_wgs_files()
genome_gr = EBVhelpR::load_wgs_reference_genome()
viral_enrich_df = EBVhelpR::load_wgs_count_summary(wgs_files_df, genome_gr)

wgs_prof_df = readRDS("../../projects_ashley/EBV_DLBCL/paper1_enhanced_EBV_detection/despiked_wgs_profiles.Rds")
wgs_prof_df$sample = factor(wgs_prof_df$sample, levels = final_ids)

ggplot(wgs_prof_df, aes(x = x, y = sample, fill = log10(y))) +
    geom_raster()

ggplot(wgs_prof_df, aes(x = x, y = sample, fill = y_norm)) +
    geom_raster()


wgs_avg_df = readRDS("../../projects_ashley/EBV_DLBCL/paper1_enhanced_EBV_detection/despiked_wgs_average_profiles.Rds")
ggplot(wgs_avg_df, aes(x = y+1, y = sample, fill = EBER_status)) +
    geom_col() +
    scale_fill_manual(values = get_colors_EBER_status()) +
    scale_x_log10()


ggplot(wgs_avg_df, aes(x = y+1, y = sample_reorder, fill = EBER_status)) +
    geom_col() +
    scale_fill_manual(values = get_colors_EBER_status()) +
    scale_x_log10()
# final_ids = rev(final_ids)

res_4plex$CTEBV_1$counts

plot_counts_barplot = function(res){
    probe_counts.wgs = .prep_probe_counts(res)

    ggplot(probe_counts.wgs, aes(x = 100*fraction, y = sample_id, fill = EBER_status)) + geom_col() +
        facet_wrap(~probe, nrow = 1) +
        scale_fill_manual(values = EBVhelpR::get_colors_EBER_status()) +
        labs(y = "Highest Percent Positive")
}

.prep_probe_counts = function(res){
    counts = lapply(res, function(x){
        x$counts
    }) %>% bind_rows(.id = "sample_id")
    counts = counts %>% pivot_wider(id_cols = c("sample_id", "unique_id"), names_from = "type", values_from = "N", values_fill = 0)


    meta_df = CellQuery()@meta_data_df %>% tibble
    counts = merge(counts,  meta_df)

    counts = counts %>% mutate(fraction = any / all) %>% arrange(fraction)
    counts$sample_id = factor(counts$sample_id, levels = counts$sample_id)

    ggplot(counts, aes(x = any/all, y = sample_id, fill = EBER_status)) + geom_col()

    counts.wgs = counts %>% filter(sample_id %in% final_ids)
    counts.wgs$sample_id = factor(counts.wgs$sample_id, levels = final_ids)

    #TODO


    ggplot(counts.wgs, aes(x = any/all, y = sample_id, fill = EBER_status)) + geom_col()


    piv_cols = c("sample_id", "unique_id", "all", "EBER_status", "sample_type")
    probe_counts.wgs = counts.wgs %>% select(-fraction) %>% pivot_longer(cols = -all_of(piv_cols), names_to = "probe", values_to = "count")
    probe_counts.wgs = probe_counts.wgs %>% mutate(fraction = count / all)

    probe_counts.wgs = probe_counts.wgs %>% filter(probe != "any")
    probe_counts.wgs
}

plot_counts_boxplot = function(res){
    probe_counts.wgs = .prep_probe_counts(res)

    ggplot(probe_counts.wgs %>% filter(probe != "any"), aes(x = EBER_status, y = fraction*100, fill = EBER_status)) + geom_boxplot() +
        facet_wrap(~probe, nrow = 1) +
        scale_fill_manual(values = EBVhelpR::get_colors_EBER_status()) +
        labs(y = "Highest Percent Positive")



}

res = res_4plex

plot_corr_plot = function(res){
    probe_counts.wgs = .prep_probe_counts(res)
    v_df = viral_enrich_df %>% select(sample_id, viral_enrichment)
    probe_counts.wgs = merge(probe_counts.wgs, v_df)
    ggplot(probe_counts.wgs, aes(x = viral_enrichment, y = fraction*100, color = EBER_status, group = 1)) +
        geom_point() +
        scale_x_log10() +
        scale_y_log10() +
        facet_wrap(~probe) +
        geom_smooth(method = "lm", formula = y ~ x) +
        ggpmisc::stat_poly_eq(aes(label = after_stat(rr.label)),
                              formula = y ~ x,
                              parse = TRUE)
}

probe_counts.wgs = .prep_probe_counts(res)
probe_counts.wgs$metric = "best_view"
cq = CellQuery(EBV_ASSAY_TYPES$RNAScope_4plex)
sum_df = cq@summary_df
sum_df$metric = "halo"

sum_df %>% head
sum_df = sum_df %>% select(sample_id, probe = combo, Positive_Percent, metric, EBER_status) %>% filter(!grepl("_", probe))
probe_counts.wgs %>% head
probe_counts.wgs = probe_counts.wgs %>% select(sample_id, EBER_status, probe, fraction, metric) %>% mutate(Positive_Percent = fraction * 100)
probe_counts.wgs$fraction = NULL

comb_df = rbind(sum_df, probe_counts.wgs)
comb_df = comb_df %>% filter(sample_id %in% final_ids)
comb_df$sample_id = factor(comb_df$sample_id, levels = final_ids)
comb_df$metric = factor(comb_df$metric, levels = c("halo", "best_view"))
ggplot(comb_df, aes(x = Positive_Percent, y = sample_id, fill = EBER_status)) +
    geom_col() +
    scale_fill_manual(values = EBVhelpR::get_colors_EBER_status()) +
    facet_wrap(~metric)
ggsave("best_barplot_combo.png", width = 13, height = 8.5)

p_bar = plot_counts_barplot(res_4plex)
ggsave(plot = p_bar, "best_barplot.png", width = 9, height = 7.5)
p_box = plot_counts_boxplot(res_4plex)
plot_corr_plot(res_4plex)

theme_set(ggpubr::theme_pubr())

plot_counts(res_4plex)
plot_counts(res_4plex) + coord_cartesian(xlim = c(0, .2))
# plot_counts(res_P)

counts_4plex = lapply(res_4plex, function(x){
    x$counts
}) %>% bind_rows(.id = "sample_id")
counts_4plex = counts_4plex %>% pivot_wider(id_cols = c("sample_id", "unique_id"), names_from = "type", values_from = "N", values_fill = 0)


meta_df = CellQuery()@meta_data_df %>% tibble
counts_4plex = merge(counts_4plex,  meta_df)

counts_4plex = counts_4plex %>% mutate(fraction = any / all) %>% arrange(fraction)
counts_4plex$sample_id = factor(counts_4plex$sample_id, levels = counts_4plex$sample_id)

ggplot(counts_4plex, aes(x = any/all, y = sample_id, fill = EBER_status)) + geom_col()

counts_4plex.wgs = counts_4plex %>% filter(sample_id %in% final_ids)
counts_4plex.wgs$sample_id = factor(counts_4plex.wgs$sample_id, levels = final_ids)

ggplot(counts_4plex.wgs, aes(x = any/all, y = sample_id, fill = EBER_status)) + geom_col()

#### view at best ####
sel_id = "D_EB_16"
sel_id = "D_EB_18"

assay_type = EBV_ASSAY_TYPES$RNAScope_4plex
plot_best_view = function(sel_id, assay_type){
    cq = CellQuery(assay_type = assay_type)
    cq = set_selected_unique_ids(cq, sel_id)
    cds = load_query_cell_data(cq)

    cq@tiff_paths_df
    tiff_df = get_query_tiff_paths_df(cq)
    xrng = res_4plex[[sel_id]]$hits$x %>% range
    yrng = res_4plex[[sel_id]]$hits$y %>% range

    hit_df = res_4plex[[sel_id]]$hits

    res_4plex[[sel_id]]$center
    view_rect = TiffRect(min(xrng), xmax = max(xrng), ymin = min(yrng), ymax = max(yrng))
    chan_names = EBV_CHANNELS[[cq@assay_type]]
    p_slide = tryCatch({
        tiff_dat = fetchTiffData(tiff_df$tiff_file[1], channel_names = chan_names)
        p_slide = tiff_dat@plots$normalized
        p_slide = p_slide %>% rect_annotate(view_rect, color = "green")
        p_slide
    }, error = function(e){
        ggplot() + labs(title = "no data")
    })

    # precalc_max = tiff_dat@precalc_max
    # precalc_max = precalc_max %>% mutate(max_value = ifelse(channel == 2, 20, max_value))
    # precalc_max = precalc_max %>% mutate(max_value = ifelse(channel == 6, 10, max_value))
    tiff_data_wide = fetchTiffData(tiff_df$tiff_file[1], rect = view_rect, channel_names = chan_names)
    tiff_data_wide.rgb = fetchTiffData.rgb(tiff_df$tiff_file[1], rect = view_rect,
                                           channel_names = chan_names, green_channel = 2, blue_channel = 1, red_channel = 6)

    # p_tiff_data = fetchTiffData(
    #     tiff_df$tiff_file[1],
    #     rect = view_rect %>% rect_resize_mult(.2),
    #     channel_names = chan_names, quantile_norm = .995)
    #
    # p_tiff_data@precalc_max$max_value[6] = 120

    close_rect = view_rect %>% rect_resize_abs(300, 300)

    p_tiff_data = fetchTiffData(
        tiff_df$tiff_file[1],
        rect = close_rect,
        channel_names = chan_names, quantile_norm = .995
    )


    chan_lev = p_tiff_data@data$channel %>% levels
    p_tiff = p_tiff_data@plots$normalized


    hit_df = res_4plex[[sel_id]]$hits
    hit_df = hit_df %>% filter(!type %in% c('any'))

    hit_df$i = hit_df$x
    hit_df$j = hit_df$y
    hit_df$channel = hit_df$type
    hit_df$norm_value = 0
    hit_df$channel %>% table
    hit_df = hit_df %>% mutate(channel = ifelse(channel == "EBER", "EBER1", channel))
    hit_df = hit_df %>% mutate(channel = ifelse(channel == "all", "DAPI", channel))
    hit_df$channel = factor(hit_df$channel, levels = chan_lev)



    p_tiff + geom_point(data = hit_df, color = 'green', shape = 1, size = 2)

    hit_rects = cds %>% filter(ObjectId %in% hit_df$name)
    decode = EBV_OPAL_DECODE[[cq@assay_type]]
    decode[decode == "EBER"] = "EBER1"
    names(decode) = paste0(names(decode), "Classification")
    hit_rects = hit_rects %>% select(ObjectId, XMin, XMax, YMin, YMax, contains("Class")) %>%
        pivot_longer(cols = contains("Class"))
    hit_rects$name = factor(hit_rects$name)
    levels(hit_rects$name) = decode[levels(hit_rects$name)]
    hit_rects$i = 1
    hit_rects$j = 1
    hit_rects$norm_value = 1

    # add all rectangles in DAPI
    dapi_rects = hit_rects %>% filter(name == hit_rects$name[1])
    dapi_rects$name = "DAPI"

    hit_rects = hit_rects %>% filter(value == 1)
    # hit_rects = rbind(hit_rects, dapi_rects)

    stopifnot(all(hit_rects$name %in% chan_names))
    hit_rects$channel = hit_rects$name %>% factor(levels = chan_names)


    bg_rect_col = "#0000FF77"
    fg_rect_col = '#00FF00AA'

    p_facet_wide = tiff_data_wide@plots$normalized

    p_facet_wide_anno = p_facet_wide +
        annotate("rect", xmin = dapi_rects$XMin, xmax = dapi_rects$XMax, ymin = dapi_rects$YMin, ymax = dapi_rects$YMax, color = bg_rect_col, fill = NA, alpha = .3) +
        geom_rect(data = hit_rects, aes(xmin = XMin, xmax = XMax, ymin = YMin, ymax = YMax), color = fg_rect_col, fill = NA)

    p_rgb_wide = tiff_data_wide.rgb@plots$rgb

    p_rgb_wide_anno = p_rgb_wide +
        annotate("rect", xmin = dapi_rects$XMin, xmax = dapi_rects$XMax, ymin = dapi_rects$YMin, ymax = dapi_rects$YMax, color = bg_rect_col, fill = NA, alpha = .3) +
        geom_rect(data = hit_rects, aes(xmin = XMin, xmax = XMax, ymin = YMin, ymax = YMax), color = fg_rect_col, fill = NA)

    p_facet_close = p_tiff
    p_facet_close_anno = p_facet_close +
        annotate("rect", xmin = dapi_rects$XMin, xmax = dapi_rects$XMax, ymin = dapi_rects$YMin, ymax = dapi_rects$YMax, color = bg_rect_col, fill = NA, alpha = .3) +
        geom_rect(data = hit_rects, aes(xmin = XMin, xmax = XMax, ymin = YMin, ymax = YMax), color = fg_rect_col, fill = NA)

    tiff_dat.rgb = fetchTiffData.rgb(tiff_df$tiff_file[1], rect = close_rect,
                                     channel_names = chan_names, green_channel = 2, blue_channel = 1, red_channel = 6)

    p_rgb_close = tiff_dat.rgb@plots$rgb
    p_rgb_close_anno = p_rgb_close +
        annotate("rect", xmin = dapi_rects$XMin, xmax = dapi_rects$XMax, ymin = dapi_rects$YMin, ymax = dapi_rects$YMax, color = bg_rect_col, fill = NA, alpha = .3) +
        geom_rect(data = hit_rects %>% filter(channel == "EBER1"), aes(xmin = XMin, xmax = XMax, ymin = YMin, ymax = YMax), color = fg_rect_col, fill = NA)

    all_plots = list(
        whole_slide = p_slide,
        facet_wide = p_facet_wide,
        facet_wide_anno = p_facet_wide_anno,
        rgb_wide = p_rgb_wide,
        rgb_wide_anno = p_rgb_wide_anno,
        facet_close = p_facet_close,
        facet_close_anno = p_facet_close_anno,
        rgb_close = p_rgb_close,
        rgb_close_anno = p_rgb_close_anno
    )
    all_plots = lapply(all_plots, function(p){
        p + theme(panel.background = element_blank()) +
            labs(x = "pixel", y = "pixel")
    })
    all_plots
}

cq = CellQuery(EBV_ASSAY_TYPES$RNAScope_4plex)
tiff_df = get_query_tiff_paths_df(cq)

img_dir = "output_cell_views_05142026_v1"
found_files = dir(img_dir, pattern = "rgb_close_anno.png")
found_samples = sub("_rgb.+", "", found_files)
todo = tiff_df$sample_id
todo = setdiff(todo, found_samples)

#these are cases that flipped
todo = c(
    "D_EB_17",
    "D_EB_18",
    "D_EB_21",
    "D_EB_23",
    "D_EB_29"
)



for(sel_id in todo){
    #rgb
    w1 = 8
    h1 = 9
    #facet
    w2 = 13
    h2 = 8.5

    dir.create(img_dir, recursive = TRUE, showWarnings = FALSE)
    cache_file = file.path(img_dir, paste0(sel_id, "_data.Rds"))
    all_plots = tryCatch({
        if(file.exists(cache_file)){
            all_plots = readRDS(cache_file)
        }else{
            all_plots = plot_best_view(sel_id, assay_type = EBV_ASSAY_TYPES$RNAScope_4plex)
            saveRDS(all_plots, file.path(img_dir, paste0(sel_id, "_data.Rds")))
        }
        all_plots
    }, error = function(e){
        NULL
    })
    if(is.null(all_plots)) next
    for(name in names(all_plots)){
        w = w2
        h = h2
        if(grepl("rgb", name)){
            w = w1
            h = h1
        }
        img_file = sub("data.Rds", paste0(name, ".png"), cache_file)
        message(img_file)
        ggsave(plot = all_plots[[name]], img_file, width = w, height = h)
    }

}
plot_best_view("D_EB_7", assay_type = EBV_ASSAY_TYPES$RNAScope_4plex)


tiff_dat@plots$normalized + geom_point(data = hit_df, color = 'green', shape = 1, size = 2)
tiff_dat2 = fetchTiffData(tiff_df$tiff_file[1], channel_names = chan_names, rect = TiffRect(xmin = 5.5e4, xmax = 6.5e4, ymax = 4.5e4, ymin = 3.5e4))
tiff_dat2@plots$normalized + geom_point(data = hit_df, color = 'green', shape = 1, size = 2)

x_rng
y_rng
cds
tr = cds %>% as.TiffRect()
tr %>% TiffPlotR::rect_centers() %>% class
