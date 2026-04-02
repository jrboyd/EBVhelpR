library(EBVhelpR)
library(ggplot2)

# Example WGS workflow using EBVhelpR helper functions.
# Update `wgs_root_dir` to your local WGS analysis directory if needed.
# Set to NULL to use package defaults.
wgs_root_dir <- NULL

meta_df <- load_meta_data()

# 1) Build WGS file index and reference ranges.
wgs_files_df <- setup_wgs_files(
  wgs_root_dir = wgs_root_dir,
  meta_df = meta_df
)

genome_gr <- load_wgs_reference_genome()

# 2) Compute host/viral read summary metrics.
wgs_count_summary <- load_wgs_count_summary(
    wgs_files_df = wgs_files_df,
    genome_gr = genome_gr,
    viral_seqname = "NC_007605.1"
)

# 3) Load smoothed bigwig pileups over the EBV genome.
pileup_df <- load_wgs_bigwig_pileup(
    wgs_files_df = wgs_files_df,
    genome_gr = genome_gr,
    viral_seqname = "NC_007605.1",
    smooth_n = 50,
    mc_cores = 4
)

# 4) Make core WGS plots.
p_enrichment <- plot_wgs_viral_enrichment(wgs_count_summary)
p_reads <- plot_wgs_viral_reads(wgs_count_summary)
p_lines <- plot_wgs_pileup_lines(pileup_df, ylim = c(0, 2000))
p_lines_zoom <- plot_wgs_pileup_lines(pileup_df, ylim = c(0, 20))
p_heat <- plot_wgs_pileup_heatmap(
    pileup_df = pileup_df,
    wgs_count_summary = wgs_count_summary,
    normalize_by_sample = FALSE
)
p_heat_norm <- plot_wgs_pileup_heatmap(
    pileup_df = pileup_df,
    wgs_count_summary = wgs_count_summary,
    normalize_by_sample = TRUE
)

colors_EBER_status = seqsetvis::safeBrew(wgs_count_summary$EBER_status)
colors_EBER_status = as.list(colors_EBER_status)
colors_EBER_status$Negative = "cornflowerblue"
colors_EBER_status$`Not performed` = "palegreen"
colors_EBER_status$Positive = "coral"
saveRDS(colors_EBER_status, "inst/extdata/colors_EBER_status.Rds")

head(wgs_count_summary)

theme_set(ggpubr::theme_pubr())

ggplot(wgs_count_summary, aes(y = sample_id, viral_read_fraction, fill = EBER_status, color = EBER_status)) +
    geom_col() +
    scale_fill_manual(values = colors_EBER_status) +
    scale_color_manual(values = colors_EBER_status) +
    theme(axis.text.y = element_text(hjust = 0, size = 8))

sel_df = wgs_files_df %>% dplyr::filter(sample_id == "D_EB_25")
f = sel_df$bigwig_file
full_df = rtracklayer::import.bw(f)
full_df = full_df %>% as.data.frame
full_df %>% dplyr::group_by(seqnames) %>% dplyr::summarise(score = max(score)) %>% dplyr::arrange(-score)

full_df %>% dplyr::filter(seqnames == "chr2") %>% dplyr::mutate(x = (start + end)/2)
plot_df = full_df %>% dplyr::filter(seqnames == "chr2")
ggplot(plot_df, aes(x = start, xend = end, y = score, yend = score)) + geom_segment()
ggplot(plot_df, aes(x = start, y = score)) + geom_path()

# Print plots in sequence for interactive dev sessions.
print(p_enrichment)
print(p_reads)
print(p_lines)
print(p_lines_zoom)
print(p_heat)
print(p_heat_norm)



# Optional save examples:
# ggsave("wgs_viral_enrichment.png", p_enrichment, width = 6, height = 8)
# ggsave("wgs_pileup_heatmap_norm.png", p_heat_norm, width = 6, height = 8)
load_rnascope_summary_files
EBV_ASSAY_TYPES$RNAScope_4plex
rscope_df = EBVhelpR::load_rnascope_summary_files()
rscope_df %>% dplyr::filter(grepl("CTEBV_15", sample_id))

rscope_df$sample_id
rscope_df$EBER_status %>% table
debug(load_phenocycler_summary_files)
phenO_df = EBVhelpR::load_phenocycler_summary_files()
phenO_df %>% subset(grepl("Neg", source)) %>% dplyr::select(EBER_status)

EBVhelpR::get_query_summary_df()

seqsetvis::ssvFeatureVenn(list(WGS = wgs_count_summary$sample_id, RSCOPE_4 = rscope_df[rscope_df$assay == EBV_ASSAY_TYPES$RNAScope_4plex,]$sample_id))
seqsetvis::ssvFeatureVenn(list(WGS = wgs_count_summary$sample_id, RSCOPE_3IF = rscope_df[rscope_df$assay == EBV_ASSAY_TYPES$`RNAScope_3plex+IF`,]$sample_id))
seqsetvis::ssvFeatureVenn(list(WGS = wgs_count_summary$sample_id, PHENO = phenO_df$sample_id))

setdiff(wgs_count_summary$sample_id, rscope_df[rscope_df$assay == EBV_ASSAY_TYPES$RNAScope_4plex,]$sample_id)
setdiff(wgs_count_summary$sample_id, rscope_df[rscope_df$assay == EBV_ASSAY_TYPES$`RNAScope_3plex+IF`,]$sample_id)
setdiff(wgs_count_summary$sample_id, phenO_df$sample_id)

setdiff(rscope_df[rscope_df$assay == EBV_ASSAY_TYPES$RNAScope_4plex,]$sample_id, wgs_count_summary$sample_id)
setdiff(rscope_df[rscope_df$assay == EBV_ASSAY_TYPES$`RNAScope_3plex+IF`,]$sample_id, wgs_count_summary$sample_id)
setdiff(phenO_df$sample_id, wgs_count_summary$sample_id)

phenO_df %>% dplyr::filter(is.na(sample_id))
