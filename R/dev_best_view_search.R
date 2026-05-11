library(EBVhelpR)
library(TiffPlotR)
library(tidyverse)

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

load_cell_source_files()
cq = EBVhelpR::CellQuery(EBV_ASSAY_TYPES$RNAScope_4plex)
todo_ids = cq@selected_unique_ids



count_in_best = function(sel_id){

    cq_sel = EBVhelpR::set_selected_unique_ids(cq, sel_id)
    cds = load_query_cell_data(cq_sel)
    message(sel_id, " found ", nrow(cds), " cells")


    tiff_df = EBVhelpR::get_query_tiff_paths_df(cq_sel)

    decode_opal_colnames = function(cq.i, cell_df = NULL){
        if(is.null(cell_df)){
            cell_df = load_query_cell_data(cq.i)
        }
        opal_codes = EBV_OPAL_DECODE[[cq.i@assay_type]]
        for(name in names(opal_codes)){
            colnames(cell_df) = sub(name, paste0(opal_codes[name], "_"), colnames(cell_df))
        }
        colnames(cell_df)
        cell_df
    }

    cds = decode_opal_colnames(cq_sel, cds)
    cds = cds %>% select(unique_id, ObjectId, XMin, XMax, YMin, YMax, contains("Classification"))


    message("classifying cells...")
    centers_all = cds %>% as.TiffRect(name_col = "ObjectId") %>% rect_centers() %>% tibble()
    # debug(as.TiffRect)
    centers_any = cds %>% filter(if_any(contains("Classification"), ~ (. == 1))) %>% as.TiffRect(name_col = "ObjectId") %>% rect_centers() %>% tibble()
    class_cns = cds %>% select(contains("Classification")) %>% colnames
    names(class_cns) = sub("_Class.+", "", class_cns)
    cn = class_cns[1]
    centers_each = lapply(class_cns, function(cn){
        hit = cds %>% filter(!!sym(cn) == 1)
        if(nrow(hit) == 0) return(NULL)
        # hit %>% as.TiffRect(name_col = "ObjectId") %>% rect_centers() %>% tibble()
        hit %>% mutate(x = (XMin + XMax)/2, y = (YMin + YMax)/2, name = ObjectId, .keep = "none")
    })
    centers_todo = c(centers_each, list(all = centers_all, any = centers_any))
    drop = sapply(centers_todo, is.null)
    centers_todo = centers_todo[!drop]

    block_size = 5e3

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
    target_cell_count = 5e3
    dist_df = centers_todo$all %>% mutate(distance = sqrt((center_point$x - x)^2 + (center_point$y - y)^2))
    dist_df = dist_df %>% mutate(rnk = rank(distance, ties.method = "first")) %>% arrange(rnk)
    best_cells_df = dist_df %>% filter(rnk <= target_cell_count)
    best_cell_ids = best_cells_df$name

    centers_nearest = lapply(centers_todo, function(x){
        x %>% filter(name %in% best_cell_ids)
    })

    ch_res = chull(centers_nearest$all$x, centers_nearest$all$y)
    #area of chull?
    ch_poly = tibble(x = centers_nearest$all$x[ch_res],
                     y = centers_nearest$all$y[ch_res]
    )

    hits_df = bind_rows(centers_nearest, .id = "type")
    count_df = hits_df %>% group_by(type) %>% summarise(N = length(type))
    count_df$unique_id = sel_id
    count_df
}

count_res = lapply(todo_ids, count_in_best)

rect_test_contains(count_rects)

x_rng
y_rng
cds
tr = cds %>% as.TiffRect()
tr %>% TiffPlotR::rect_centers() %>% class
