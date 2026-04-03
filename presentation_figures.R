library(EBVhelpR)
library(ggplot2)
library(tidyverse)
# Example WGS workflow using EBVhelpR helper functions.
# Update `wgs_root_dir` to your local WGS analysis directory if needed.
# Set to NULL to use package defaults.

# 01 WGS viral enrichment and for all samples. 1 sample with no EBER status omitted, "not performed". Includes profiles along viral genome.
# 02 overlap of samples completed for WGS, RNAscope1+2 and Phenocycler. Final selection of 30 patients with complete set of data.
# 03 replot of WGS for selected final patient set
# 04 summaries for RNAscope1+2 and Phenocycler. This is the same percent positive data Kyra generated.
# 05 final select of scope data for cell viewing
# 06 correlation of scope data with WGS

res_dir = "output_presentation_04092026_v2"
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

#### group WGS initial plots ####
increase_plot_group()
print(p_enrichment)
my_saveplot(p_enrichment, "wgs_enrichment_barplot", width = 11, height = 8)
print(pg_lines)
my_saveplot(pg_lines, "wgs_profiles_lineplot", width = 11, height = 8)
print(p_heat)
my_saveplot(p_heat, "wgs_profiles_heatmap", width = 11, height = 8)
print(p_heat_norm)
my_saveplot(p_heat_norm, "wgs_profiles_heatmap_normalized", width = 11, height = 8)

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

#### group overlap ids ####
increase_plot_group()
upr = seqsetvis::ssvFeatureUpset(all_ids, return_UpSetR = TRUE)

pdf(res_file("assay_overlap_upset.pdf"), width = 6, height = 5)
print(upr)
dev.off()
png(res_file("assay_overlap_upset.png"), width = 6, height = 5, units = "in", res = 200)
print(upr)
dev.off()

membs = seqsetvis::ssvMakeMembTable(all_ids)
in_all = apply(membs, 1, all)
final_ids = rownames(membs)[in_all]
saveRDS(final_ids, res_file("final_ids.Rds"))


#### apply final ids ####
wgs_count_summary.re = wgs_count_summary %>% filter(sample_id %in% final_ids)
wgs_count_summary.re$sample_id = wgs_count_summary.re$sample_id %>% droplevels()
pileup_df.re = pileup_df %>% filter(sample %in% final_ids)
cq_r4 = set_selected_unique_ids(cq_r4, final_ids)
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

#### group replot WGS ####
increase_plot_group()
print(p_enrichment.re)
my_saveplot(p_enrichment.re, "wgs_enrichment_barplot", width = 11, height = 8)
print(pg_lines.re)
my_saveplot(pg_lines.re, "wgs_profiles_lineplot", width = 11, height = 8)
print(p_heat.re)
my_saveplot(p_heat.re, "wgs_profiles_heatmap", width = 11, height = 8)
print(p_heat_norm.re)
my_saveplot(p_heat_norm.re, "wgs_profiles_heatmap_normalized", width = 11, height = 8)

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
        scale_fill_manual(values = colors_EBER_status) +
        labs(x = "Percent Positive", y = "", fill = "EBER status")

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
        scale_fill_manual(values = colors_EBER_status) +
        labs(x = "Probe", y = "Percent Positive", fill = "EBER status")

}


p_sum_bars_r4 = cq_plot_summary_bars(cq_r4)
p_sum_bars_r4

p_sum_bars_r3i = cq_plot_summary_bars(cq_r3i)
p_sum_bars_r3i


p_sum_box_r4 = cq_plot_summary_boxplot(cq_r4)
p_sum_box_r4
p_sum_box_r3i = cq_plot_summary_boxplot(cq_r3i)
p_sum_box_r3i


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

p_sum_bars_p = ggplot(pos_p, aes(x = value, y = sample_id, fill = EBER_status)) +
    geom_col() +
    facet_wrap(~name, scales = "free_x", nrow = 2) +
    scale_fill_manual(values = colors_EBER_status) +
    scale_x_continuous(breaks = scales::pretty_breaks(2)) +
    labs(x = "Percent Positive", y = "", fill = "EBER status") +
    theme(axis.text.y = element_text(size= 8), axis.text.x = element_text(size = 8))

p_sum_bars_p


p_sum_box_p = ggplot(pos_p, aes(x = name, y = value, fill = EBER_status)) +
    geom_boxplot() +
    # facet_wrap(~name, scales = "free_x", nrow = 1) +
    scale_fill_manual(values = colors_EBER_status) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(y = "Percent Positive", x = "Probe", fill = "EBER status")
p_sum_box_p

#### group scope summaries ####

increase_plot_group()
my_saveplot(p_sum_bars_r4, name = "RNAScope_4plex_barplot", width = 11, height = 7.2)
my_saveplot(p_sum_bars_r3i, name = "RNAScope_3plexIF_barplot", width = 11, height = 7.2)
my_saveplot(p_sum_bars_p, name = "Phenocycler_barplot", width = 11, height = 8.3)

my_saveplot(p_sum_box_r4, name = "RNAScope_4plex_boxplot", width = 8.4, height = 5)
my_saveplot(p_sum_box_r3i, name = "RNAScope_3plexIF_boxplot", width = 8.4, height = 5)
my_saveplot(p_sum_box_p, name = "Phenocycler_boxplot", width = 7.5, height = 6.4)


#### cell images ####
increase_plot_group()
saveRDS(cq_r4, res_file("cq_r4.Rds"))
saveRDS(cq_r3i, res_file("cq_r3i.Rds"))
saveRDS(cq_p, res_file("cq_p.Rds"))

qsum = get_query_summary_df(cq_r4)
qsum %>% head

xy_df = merge(qsum, wgs_count_summary.re %>% select(sample_id, viral_enrichment), by = "sample_id")
xy_df$sample_id = factor(xy_df$sample_id, levels = id_lev)
xy_df %>% filter(!grepl("Unstai", combo))

# install.packages("ggpmisc")
library(ggpmisc)
increase_plot_group()
p_perc_vs_enrich = ggplot(xy_df, aes(x = viral_enrichment, y = Positive_Percent, color = EBER_status)) +
    geom_point() +
    facet_wrap(~combo, scales = "free_y") +
    scale_x_log10() +
    scale_y_log10() +
    stat_poly_line(method = "lm", formula = y ~ x) +
    stat_poly_eq(aes(label = ..rr.label..),
                 formula = y ~ x, parse = TRUE, label.y = c(.8, .95))
p_perc_vs_enrich
my_saveplot(p_perc_vs_enrich, name = "RNAscope4plex_vs_WGS_grouped", width = 9.3, height = 7.5)

p_perc_vs_enrich.grouped = ggplot(xy_df, aes(x = viral_enrichment, y = Positive_Percent, color = EBER_status, group = combo)) +
    geom_point() +
    facet_wrap(~combo, scales = "free_y") +
    scale_x_log10() +
    scale_y_log10() +
    stat_poly_line(method = "lm", formula = y ~ x) +
    stat_poly_eq(aes(label = ..rr.label..),
                 formula = y ~ x, parse = TRUE, label.y = c(.8, .95))
my_saveplot(p_perc_vs_enrich.grouped, name = "RNAscope4plex_vs_WGS_ungrouped", width = 9.3, height = 7.5)

# 3plexIF
qsum = get_query_summary_df(cq_r3i)
qsum %>% head

xy_df = merge(qsum, wgs_count_summary.re %>% select(sample_id, viral_enrichment), by = "sample_id")
xy_df$sample_id = factor(xy_df$sample_id, levels = id_lev)
xy_df %>% filter(!grepl("Unstai", combo))

# install.packages("ggpmisc")
# library(ggpmisc)
# increase_plot_group()
p_perc_vs_enrich = ggplot(xy_df, aes(x = viral_enrichment, y = Positive_Percent, color = EBER_status)) +
    geom_point() +
    facet_wrap(~combo, scales = "free_y") +
    scale_x_log10() +
    scale_y_log10() +
    stat_poly_line(method = "lm", formula = y ~ x) +
    stat_poly_eq(aes(label = ..rr.label..),
                 formula = y ~ x, parse = TRUE, label.y = c(.8, .95))
p_perc_vs_enrich
my_saveplot(p_perc_vs_enrich, name = "RNAscope3plexIF_vs_WGS_grouped", width = 9.3, height = 7.5)

p_perc_vs_enrich.grouped = ggplot(xy_df, aes(x = viral_enrichment, y = Positive_Percent, color = EBER_status, group = combo)) +
    geom_point() +
    facet_wrap(~combo, scales = "free_y") +
    scale_x_log10() +
    scale_y_log10() +
    stat_poly_line(method = "lm", formula = y ~ x) +
    stat_poly_eq(aes(label = ..rr.label..),
                 formula = y ~ x, parse = TRUE, label.y = c(.8, .95))
my_saveplot(p_perc_vs_enrich.grouped, name = "RNAscope3plexIF_vs_WGS_ungrouped", width = 9.3, height = 7.5)

#Pheno
# 3plexIF
qsum = get_query_summary_df(cq_p)
pos_p = qsum_p %>% select(sample_id, contains("Cells") & contains("%")) %>%
    pivot_longer(cols = !sample_id)
pos_p = pos_p %>% mutate(name = sub("% ", "", name)) %>% mutate(name = sub(" Positive Cells", "", name))
pos_p = merge(meta_df, pos_p)
pos_p$sample_id = factor(pos_p$sample_id, levels = id_lev)
head(pos_p)
qsum =pos_p %>% rename(Positive_Percent = value, combo = name)

xy_df = merge(qsum, wgs_count_summary.re %>% select(sample_id, viral_enrichment), by = "sample_id")
xy_df$sample_id = factor(xy_df$sample_id, levels = id_lev)
xy_df %>% filter(!grepl("Unstai", combo))


pheno_combos = xy_df$combo %>% unique
n_per = 9
todo = split(pheno_combos, ceiling(seq_along(pheno_combos)/n_per))
xy_df.full = xy_df
i = 1
for(i in seq_along(todo)){
    sel_combos = todo[[i]]
    xy_df = xy_df.full %>% filter(combo %in% sel_combos)
    n_combos = xy_df$combo %>% unique %>% length

    if(n_combos < 9){
        # add placeholders so plot sizes are consistent
        n_add = 9 - n_combos
        xy_df$combo = factor(xy_df$combo, levels = c(unique(xy_df$combo), paste("empty", LETTERS[seq(n_add)])))
    }

    # install.packages("ggpmisc")
    # library(ggpmisc)
    # increase_plot_group()
    p_perc_vs_enrich = ggplot(xy_df, aes(x = viral_enrichment, y = Positive_Percent, color = EBER_status)) +
        geom_point() +
        facet_wrap(~combo, scales = "free_y", ncol = 3, nrow = 3, drop = FALSE) +
        scale_x_log10() +
        scale_y_log10() +
        stat_poly_line(method = "lm", formula = y ~ x) +
        stat_poly_eq(aes(label = ..rr.label..),
                     formula = y ~ x, parse = TRUE, label.y = c(.8, .95))
    p_perc_vs_enrich
    f_suffix = paste0(".", i, "_of_", length(todo))
    my_saveplot(p_perc_vs_enrich, name = paste0("Phenocycler_vs_WGS_grouped", f_suffix), width = 9.3, height = 7.5)

    p_perc_vs_enrich.grouped = ggplot(xy_df, aes(x = viral_enrichment, y = Positive_Percent, color = EBER_status, group = combo)) +
        geom_point() +
        facet_wrap(~combo, scales = "free_y", ncol = 3, nrow = 3, drop = FALSE) +
        scale_x_log10() +
        scale_y_log10() +
        stat_poly_line(method = "lm", formula = y ~ x) +
        stat_poly_eq(aes(label = ..rr.label..),
                     formula = y ~ x, parse = TRUE, label.y = c(.8, .95))
    my_saveplot(p_perc_vs_enrich.grouped, name = paste0("Phenocycler_vs_WGS_ungrouped", f_suffix), width = 9.3, height = 7.5)
}
