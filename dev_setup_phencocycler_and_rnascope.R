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

query_df = cell_df %>% filter(
    Opal520Classification == 1
)
query_df.l = split(query_df, query_df$sample_id)
n_cells = 9
query_df.l = lapply(query_df.l, function(x){
    x[sample(nrow(x), size = n_cells),]
})

library(TiffPlotR)
rects_from_df = function(df){
    lapply(seq(nrow(df)), function(i){
        TiffRect(df$XMin[i], df$XMax[i], df$YMin[i], df$YMax[i])
    })

}
rects.l = lapply(query_df.l, rects_from_df)
head(query_df.l$CTEBV_11)

query_df.l$CTEBV_11$XMin

# tiff image paths
image_dir = "Z:/FUSION DATA/AshleyVolaric/RNAScopeRound1/"
image_files = sapply(names(query_df.l), function(sample_id){
    dir(image_dir, pattern = gsub("_", "", sample_id), full.names = TRUE)
}) %>% as.list

img_f = image_files$CTEBV_11
rects = rects.l$CTEBV_11

r = rects[[1]]
r_fetch = r %>% rect_resize_mult(2)
img_res = fetchTiffData(img_f, r_fetch)
# undebug(fetchTiffData.rgb)
img_res = fetchTiffData.rgb(img_f, r_fetch)
img_res@plots$normalized
rect_annotate(img_res@plots$normalized, r)
img_res@plots$normalized %>% rect_annotate(., r)

#
# p = img_res@plots$normalized
# rect = r
# rect_df <- data.frame(xmin = rect@xmin, xmax = rect@xmax, ymin = rect@ymin, ymax = rect@ymax)
# p + ggplot2::geom_rect(data = rect_df, mapping = ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE, color = "red")
