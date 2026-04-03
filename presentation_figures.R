library(EBVhelpR)
library(ggplot2)
library(tidyverse)
# Example WGS workflow using EBVhelpR helper functions.
# Update `wgs_root_dir` to your local WGS analysis directory if needed.
# Set to NULL to use package defaults.

res_dir = "output_presentation_04092026_v1"
res_file = function(f){
  out_f = file.path(res_dir, f)
  dir.create(dirname(out_f), recursive = TRUE, showWarnings = FALSE)
  out_f
}
my_saveplot = function(plot, name, width, height){
  ggsave(res_file(paste0(name, ".png")), plot, width = width, height = height)
  ggsave(res_file(paste0(name, ".pdf")), plot, width = width, height = height)
}

#### WGS setup ####

wgs_root_dir <- NULL

meta_df <- load_meta_data()
colors_EBER_status = get_colors_EBER_status()
theme_set(ggpubr::theme_pubr())

# 1) Build WGS file index and reference ranges.
wgs_files_df <- setup_wgs_files(
  wgs_root_dir = wgs_root_dir,
  meta_df = meta_df
)

#remove Not Performed sample
wgs_files_df = wgs_files_df %>% filter(EBER_status %in% c("Positive", "Negative"))

genome_gr <- load_wgs_reference_genome()

# 2) Compute host/viral read summary metrics.
wgs_count_summary <- load_wgs_count_summary(
  wgs_files_df = wgs_files_df,
  genome_gr = genome_gr,
  viral_seqname = "NC_007605.1"
)

# 3) Load smoothed bigwig pileups over the EBV genome.
cache_dir = "presentation_cache"
pileup_cache = file.path(cache_dir, "pileup_df.Rds")

if(file.exists(pileup_cache)){
  pileup_df = readRDS(pileup_cache)
}else{
  pileup_df <- load_wgs_bigwig_pileup(
    wgs_files_df = wgs_files_df,
    genome_gr = genome_gr,
    viral_seqname = "NC_007605.1",
    smooth_n = 50,
    mc_cores = 4
  )
  dir.create(dirname(pileup_cache))
  saveRDS(pileup_df, pileup_cache)
}

#### WGS make plots ####

# 4) Make core WGS plots.
p_enrichment <- plot_wgs_viral_enrichment(wgs_count_summary) +
  scale_fill_manual(values = colors_EBER_status) +
  scale_color_manual(values = colors_EBER_status) +
  theme(axis.text.y = element_text(hjust = 0, size = 8))
# p_reads <- plot_wgs_viral_reads(wgs_count_summary)
p_lines <- plot_wgs_pileup_lines(pileup_df, ylim = c(0, 2000))
p_lines_zoom <- plot_wgs_pileup_lines(pileup_df, ylim = c(0, 20))

pg_lines = cowplot::plot_grid(p_lines, p_lines_zoom, nrow = 1)

p_heat <- plot_wgs_pileup_heatmap(
  pileup_df = pileup_df,
  wgs_count_summary = wgs_count_summary,
  normalize_by_sample = FALSE,
  status_colors = colors_EBER_status
) +
  scale_fill_viridis_c() +
  theme(axis.text.y = element_text(hjust = 0, size = 8))
p_heat
p_heat_norm <- plot_wgs_pileup_heatmap(
  pileup_df = pileup_df,
  wgs_count_summary = wgs_count_summary,
  normalize_by_sample = TRUE,
  status_colors = colors_EBER_status
) +
  scale_fill_viridis_c() +
  theme(axis.text.y = element_text(hjust = 0, size = 8))
p_heat_norm

#### WGS save plots ####

print(p_enrichment)
my_saveplot(p_enrichment, "01_wgs_enrichment_barplot", width = 11, height = 8)
print(pg_lines)
my_saveplot(pg_lines, "01_wgs_profiles_lineplot", width = 11, height = 8)
print(p_heat)
my_saveplot(p_heat, "01_wgs_profiles_heatmap", width = 11, height = 8)
print(p_heat_norm)
my_saveplot(p_heat_norm, "01_wgs_profiles_heatmap_normalized", width = 11, height = 8)

ids_wgs = levels(wgs_count_summary$sample_id)

#### RNAscope 4plex ####
cq_r4 = CellQuery(assay_type = EBV_ASSAY_TYPES$RNAScope_4plex)
sel_r4 = cq_r4@summary_df %>% subset(sample_type == "patient" & probe_control == "" & EBER_status %in% wgs_count_summary$EBER_status)
ids_r4 = sel_r4$unique_id %>% unique()
cq_r4 = set_selected_unique_ids(cq_r4, ids_wgs)

#### RNAscope 3plexIF ####
cq_r3i = CellQuery(assay_type = EBV_ASSAY_TYPES$`RNAScope_3plex+IF`)
sel_r3i = cq_r3i@summary_df %>% subset(sample_type == "patient" & probe_control == "" & EBER_status %in% wgs_count_summary$EBER_status)
ids_r3i = sel_r3i$unique_id %>% unique()
cq_r3i = set_selected_unique_ids(cq_r3i, ids_wgs)

#### Phenocycler ####
cq_p = CellQuery(assay_type = EBV_ASSAY_TYPES$Phenocycler)
sel_p = cq_p@summary_df %>% subset(sample_type == "patient" & probe_control == "" & EBER_status %in% wgs_count_summary$EBER_status)
ids_p = sel_p$unique_id %>% unique()
cq_p = set_selected_unique_ids(cq_p, ids_wgs)


#### select final ids ####

all_ids = list(
  WGS = ids_wgs,
  `RNAscope 4plex` = ids_r4,
  `RNAscope 3plexIF` = ids_r3i,
  Phenocycler = ids_p
)

upr = seqsetvis::ssvFeatureUpset(all_ids, return_UpSetR = TRUE)
pdf(res_file("02_assay_overlap_upset.pdf"), width = 6, height = 5)
print(upr)
dev.off()
png(res_file("02_assay_overlap_upset.png"), width = 6, height = 5, units = "in", res = 200)
print(upr)
dev.off()

membs = seqsetvis::ssvMakeMembTable(all_ids)
in_all = apply(membs, 1, all)
final_ids = rownames(membs)[in_all]


#### apply final ids ####
wgs_count_summary.re = wgs_count_summary %>% filter(sample_id %in% final_ids)
wgs_count_summary.re$sample_id = wgs_count_summary.re$sample_id %>% droplevels()
pileup_df.re = pileup_df %>% filter(sample %in% final_ids)
cq_r4 = set_selected_unique_ids(cq_r4, final_ids)
get_full_cell_data(cq_r4)
cq_r4@selected_unique_ids
get_query_summary_df(cq_r4)
cq_r3i = set_selected_unique_ids(cq_r3i, final_ids)
cq_p = set_selected_unique_ids(cq_p, final_ids)

id_lev = levels(wgs_count_summary.re$sample_id)

#### replot WGS ####
p_enrichment.re <- plot_wgs_viral_enrichment(wgs_count_summary.re) +
  scale_fill_manual(values = colors_EBER_status) +
  scale_color_manual(values = colors_EBER_status) +
  theme(axis.text.y = element_text(hjust = 0, size = 8))

p_lines.re <- plot_wgs_pileup_lines(pileup_df.re, ylim = c(0, 2000))
p_lines_zoom.re <- plot_wgs_pileup_lines(pileup_df.re, ylim = c(0, 20))

pg_lines.re = cowplot::plot_grid(p_lines.re, p_lines_zoom.re, nrow = 1)

p_heat.re <- plot_wgs_pileup_heatmap(
  pileup_df = pileup_df.re,
  wgs_count_summary = wgs_count_summary.re,
  normalize_by_sample = FALSE,
  status_colors = colors_EBER_status
) +
  scale_fill_viridis_c() +
  theme(axis.text.y = element_text(hjust = 0, size = 8))
p_heat.re
p_heat_norm.re <- plot_wgs_pileup_heatmap(
  pileup_df = pileup_df.re,
  wgs_count_summary = wgs_count_summary.re,
  normalize_by_sample = TRUE,
  status_colors = colors_EBER_status
) +
  scale_fill_viridis_c() +
  theme(axis.text.y = element_text(hjust = 0, size = 8))
p_heat_norm.re

#### save replot WGS ####
print(p_enrichment.re)
my_saveplot(p_enrichment.re, "03_wgs_enrichment_barplot", width = 11, height = 8)
print(pg_lines.re)
my_saveplot(pg_lines.re, "03_wgs_profiles_lineplot", width = 11, height = 8)
print(p_heat.re)
my_saveplot(p_heat.re, "03_wgs_profiles_heatmap", width = 11, height = 8)
print(p_heat_norm.re)
my_saveplot(p_heat_norm.re, "03_wgs_profiles_heatmap_normalized", width = 11, height = 8)

#### RNAscope 4plex plots ####
cq_plot_summary_bars = function(cq){
  qsum = get_query_summary_df(cq)
  qcell_files = get_query_cell_files_df(cq)
  qtiff = get_query_tiff_paths_df(cq)
  
  qsum$sample_id = factor(qsum$sample_id, levels = id_lev)
  head(qsum)
  qsum = qsum %>% subset(!grepl("_", combo))
  ggplot(qsum, aes(x = Positive_Percent, y = sample_id, fill = EBER_status)) +
    geom_col() +
    facet_wrap(~combo, scales = "free_x", nrow = 1) +
    scale_fill_manual(values = colors_EBER_status)
  
}

cq_plot_summary_boxplot = function(cq){
  qsum = get_query_summary_df(cq)
  qcell_files = get_query_cell_files_df(cq)
  qtiff = get_query_tiff_paths_df(cq)
  
  qsum$sample_id = factor(qsum$sample_id, levels = id_lev)
  head(qsum)
  qsum = qsum %>% subset(!grepl("_", combo))
  ggplot(qsum, aes(x = combo, y = Positive_Percent, fill = EBER_status)) +
    geom_boxplot() +
    # facet_wrap(~combo, scales = "free_x", nrow = 1) +
    scale_fill_manual(values = colors_EBER_status)
  
}

cq_plot_summary_bars(cq_r4)
cq_plot_summary_bars(cq_r3i)

cq_plot_summary_boxplot(cq_r4)
cq_plot_summary_boxplot(cq_r3i)

qsum_r4 = get_query_summary_df(cq_r4)
qcell_files_r4 = get_query_cell_files_df(cq_r4)
qtiff_r4 = get_query_tiff_paths_df(cq_r4)

qsum_p = get_query_summary_df(cq_p)
pos_p = qsum_p %>% select(sample_id, contains("Cells") & contains("%")) %>%
  pivot_longer(cols = !sample_id)
pos_p = pos_p %>% mutate(name = sub("% ", "", name)) %>% mutate(name = sub(" Positive Cells", "", name))
pos_p = merge(meta_df, pos_p)
pos_p$sample_id = factor(pos_p$sample_id, levels = id_lev)
pos_p$value %>% class
ggplot(pos_p, aes(x = value, y = sample_id, fill = EBER_status)) +
  geom_col() +
  facet_wrap(~name, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = colors_EBER_status)

ggplot(pos_p, aes(x = name, y = value, fill = EBER_status)) +
  geom_boxplot() +
  # facet_wrap(~name, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = colors_EBER_status)

#### cell images ####
cq = cq_r4
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

library(TiffPlotR)
img_res = fetchTiffData.rgb(tiff_path = tiff_f, channel_names = EBV_CHANNELS[[assay_name]], blue_channel = 1, red_channel = 2, green_channel = 3)
tp = sample(nrow(cell_df), 1e3)
img_res@plots$rgb +
  annotate("point", x= cell_df$x[tp], y = cell_df$y[tp], color = "green") +
  coord_fixed()

r = TiffRect(2e4, 2.2e4, 8e3, 9e3)
img_res.z = fetchTiffData.rgb(tiff_path = tiff_f, rect = r,
                              channel_names = EBV_CHANNELS[[assay_name]], blue_channel = 1, red_channel = 2, green_channel = 3)
tp = sample(nrow(cell_df), 1e3)
img_res.z@plots$rgb +
  # annotate("point", x= cell_df$x[tp], y = cell_df$y[tp], color = "green") +
  coord_fixed(xlim = c(r@coords$xmin, r@coords$xmax), ylim = c(r@coords$ymin, r@coords$ymax))

r2 = rect_resize_mult(r, .2)
fetchTiffData.rgb(tiff_path = tiff_f, rect = r2,
                  channel_names = EBV_CHANNELS[[assay_name]], blue_channel = 1, red_channel = 4, green_channel = 5)

TiffRect()
