library(magrittr)
library(tidyverse)
library(EBVhelpR)
#TODO resolve inconsistent behaviors
# sample_id everywhere
# source vs assay confusion


EBVhelpR::write_all_package_data()


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

cq = CellQuery()
cq@summary_df

set.seed(0)
sel_ids = cq@summary_df$sample_id %>% sample(4) %>% unique

# debug(set_selected_sample_ids)
cq = set_selected_sample_ids(cq, sel_ids)

cell_df = get_query_cell_files_df(cq)
cell_df$sample_id %>% table
cell_df %>% head

cell_df = load_query_cell_data(cq)
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
