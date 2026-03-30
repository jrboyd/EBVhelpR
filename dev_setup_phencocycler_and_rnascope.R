library(magrittr)
library(tidyverse)
library(EBVhelpR)
#TODO resolve inconsistent behaviors
# sample_id everywhere
# source vs assay confusion


EBVhelpR::write_all_package_data()

undebug(load_phenocycler_summary_files)
pcycler_df = EBVhelpR::load_phenocycler_summary_files()
pcycler_df$source %>% table
pcycler_df$SampleID %>% table
pcycler_df$combo %>% table
pcycler_df$EBER_status %>% table


rscope_df = EBVhelpR::load_rnascope_summary_files()
rscope_df$source %>% table
rscope_df$SampleID %>% table
rscope_df$combo %>% table
rscope_df$EBER_status %>% table

pcycler_df$Sample
rscope_df$SampleID




theme_set(theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)))
ggplot(rscope_df, aes(x = combo, y = Positive_Percent)) +
    geom_boxplot() +
    facet_wrap(~source)

# view example cells
# debug(load_cell_source_files)
cell_info_df = EBVhelpR::load_cell_source_files()


# select samples to view
rscope_df$source %>% table
sel_df = filter(rscope_df, source == "4plex" & EBER_status == "Positive")
sel_ids = sel_df$SampleID %>% unique
# just a couple for test speed
sel_ids = sel_ids[1:2]

head(cell_info_df)
sel_cell_info_df = filter(cell_info_df, sample_id %in% sel_ids & assay == "RNAScope_4plex")
ddir = get_wrangled_cell_data_dir()
cell_files = file.path(ddir, sel_cell_info_df$file)
stopifnot(file.exists(cell_files))
names(cell_files) = sel_cell_info_df$sample_id

cell_df.l = lapply(cell_files, read.csv)
cell_df = bind_rows(cell_df.l, .id = "sample_id")
lapply(cell_df, range)

ggplot(cell_df, aes(x = sample_id, y = Opal620CellIntensity)) +
    geom_boxplot()

ggplot(cell_df, aes(x = sample_id, fill = factor(Opal620Classification))) +
    geom_bar()


cell_df$Opal620Classification %>% table
class_fraction = cell_df %>%
    select(sample_id, matches("Class")) %>%
    pivot_longer(cols = !sample_id) %>%
    group_by(sample_id, name) %>%
    summarise(N = sum(value == 1), total = length(value), fraction = N / total)

ggplot(class_fraction, aes(x = sample_id, fill = name, y = fraction)) +
    geom_col(position = "dodge")

#### look at example cells ####

library(TiffPlotR)

select_representative_cells <- function(
        cell_df,
        marker_col = "Opal520Classification",
        marker_value = 1,
        n_cells = 9,
        sample_col = "sample_id",
        seed = 1
) {
    stopifnot(marker_col %in% colnames(cell_df))
    stopifnot(sample_col %in% colnames(cell_df))

    query_df <- cell_df %>% dplyr::filter(!!dplyr::sym(marker_col) == marker_value)
    query_df.l <- split(query_df, query_df[[sample_col]])

    set.seed(seed)
    query_df.l <- lapply(query_df.l, function(x) {
        if (!nrow(x)) {
            return(x)
        }
        n_take <- min(nrow(x), n_cells)
        x[sample(nrow(x), size = n_take), ]
    })

    query_df.l
}

rects_from_df <- function(df) {
    lapply(seq(nrow(df)), function(i) {
        TiffPlotR::TiffRect(df$XMin[i], df$XMax[i], df$YMin[i], df$YMax[i])
    })
}

find_tiff_file_by_sample <- function(selected_sample_id, image_dir) {
    select_res = cq@tiff_paths_df %>% dplyr::filter(sample_id == selected_sample_id)
    if(nrow(select_res) == 0){
        stop("no matching tiff found")
    }
    if(nrow(select_res) > 1){
        stop("ambiguous tiff selection")
    }
    select_res$tiff_file[1]
}

cq = CellQuery()
cq@summary_df

set.seed(0)
sel_ids = cq@summary_df$sample_id %>% sample(4) %>% unique

# debug(set_selected_sample_ids)
cq = set_selected_sample_ids(cq, sel_ids)

cell_df = get_query_cell_files_df(cq)
cell_df$sample_id %>% table
cell_df %>% head

load_cell_files = function(cq){
    cell_info_df = get_query_cell_files_df(cq)
    to_load = cell_info_df$file
    names(to_load) = cell_info_df$sample_id

    cell_df.l = lapply(to_load, function(f){
        f = file.path(get_wrangled_cell_data_dir(), f)
        # .cell_df = read.csv(f)
        .cell_df = readr::read_csv(f)
        # .cell_df = data.table::fread(f)
        .cell_df$AlgorithmName = NULL
        .cell_df$...1 = NULL
        .cell_df
    })
    cell_df = dplyr::bind_rows(cell_df.l, .id = "sample_id")
    cell_df
}



#' Title
#'
#' @param sampled_cells
#' @param image_dir
#' @param fetch_resize_mult
#' @param max_images_per_sample
#'
#' @returns
#' @export
#'
#' @examples
fetch_representative_tiff_images <- function(
        sampled_cells,
        image_dir,
        fetch_resize_mult = 2,
        max_images_per_sample = 3
) {
    if(is.data.frame(sampled_cells)){
        stopifnot("sample_id" %in% colnames(sampled_cells))
        sampled_cells = split(sampled_cells, sampled_cells$sample_id)
    }
    sample_ids <- names(sampled_cells)

    image_files = sapply(sample_ids, find_tiff_file_by_sample, image_dir = image_dir)

    image_res <- lapply(sample_ids, function(sample_id) {
        sample_df <- sampled_cells[[sample_id]]
        img_file <- image_files[[sample_id]]

        if (is.na(img_file) || !nrow(sample_df)) {
            return(list())
        }

        rects <- rects_from_df(sample_df)
        rects <- rects[seq_len(min(length(rects), max_images_per_sample))]

        lapply(rects, function(r) {
            r_fetch <- r %>% TiffPlotR::rect_resize_mult(fetch_resize_mult)
            img_res <- TiffPlotR::fetchTiffData.rgb(img_file, r_fetch)
            img_res@plots$rgb = img_res@plots$rgb %>% TiffPlotR::rect_annotate(., r, color = "yellow")
            img_res
        })
    })
    names(image_res) <- sample_ids

    image_res
}

cell_df = load_cell_files(cq)
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

# tiff image paths
get_image_dir = function(cq){
    win_path = "Z:/FUSION DATA/AshleyVolaric"
    data_path = win_path
    if(!dir.exists(data_path)){
        stop("can't locate image directory:")
    }

    image_dir = file.path(data_path, EBVhelpR:::assay_to_project_name[cq@assay_type])
    stopifnot(file.exists(image_dir))
    image_dir

}
# image_dir = /RNAScopeRound1/"
# debug(fetch_representative_tiff_images)
image_dir = get_image_dir(cq)
tiff_sample <- fetch_representative_tiff_images(
    sampled_cells = query_df.l,
    image_dir = image_dir,
    fetch_resize_mult = 2,
    max_images_per_sample = 3
)

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
