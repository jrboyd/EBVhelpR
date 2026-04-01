library(EBVhelpR)
library(ggplot2)

# Example WGS workflow using EBVhelpR helper functions.
# Update `wgs_root_dir` to your local WGS analysis directory if needed.



# stopifnot(dir.exists(wgs_bwa_dir))
stopifnot(dir.exists(count_dir))
stopifnot(dir.exists(bigwig_dir))

meta_df <- load_meta_data()

# 1) Build reference GRanges from BAM header.
genome_gr <- load_wgs_reference_genome()

setup_wgs_files = function(){
    win_dir = "C:/Users/boydj/OneDrive - UVM Larner College of Medicine/projects_ashley/EBV_DLBCL/P2_viral_WGS"
    lin_dir = "/gpfs1/pi/avolaric/files_jrboyd/P2_viral_WGS"
    wgs_root_dir = win_dir
    if(!dir.exists(wgs_root_dir)){
        wgs_root_dir = lin_dir
        if(!dir.exists(wgs_root_dir)){
            stop("Could not locate WGS data directory.")
        }
    }

    # wgs_bwa_dir <- file.path(wgs_root_dir, "output_bwa")
    count_dir <- file.path(wgs_root_dir, "chr_read_counts")
    bigwig_dir <- file.path(wgs_root_dir, "bigwigs")
    stopifnot(dir.exists(count_dir))
    stopifnot(dir.exists(bigwig_dir))

    count_files = dir(count_dir, full.names = TRUE)
    count_df = data.frame(count_file = count_files)
    count_df = count_df %>% dplyr::mutate(sample_id = sub("_read.+", "", basename(count_file)))

    bw_files = dir(bigwig_dir, pattern = "norm", full.names = TRUE)
    bw_df = data.frame(bigwig = bw_files)
    bw_df = bw_df %>% dplyr::mutate(sample_id = sub("_dedu.+", "", basename(bigwig)))
    wgs_df = merge(count_df, bw_df, by = "sample_id", all = TRUE)
    #some samples end in extra 2, redone?
    wgs_df = wgs_df %>% dplyr::mutate(sample_id = ifelse(grepl("[0-9]_1$", sample_id), sub("_1$", "", sample_id), sample_id))
    wgs_df = wgs_df %>% dplyr::mutate(sample_id = ifelse(grepl("[0-9]_2$", sample_id), sub("_2$", "", sample_id), sample_id))
    wgs_df = merge(wgs_df, meta_df, by = "sample_id", all.x = TRUE)
    stopifnot(!any(is.na(wgs_df$EBER_status)))
    wgs_df %>% dplyr::filter(is.na(EBER_status))
    wgs_df
}

wgs_df = setup_wgs_files()

# 2) Compute host/viral read summary metrics.
wgs_count_summary <- load_wgs_count_summary(
    count_dir = count_dir,
    genome_gr = genome_gr,
    viral_seqname = "NC_007605.1",
    meta_df = meta_df
)

# 3) Load smoothed bigwig pileups over the EBV genome.
pileup_df <- load_wgs_bigwig_pileup(
    bigwig_dir = bigwig_dir,
    genome_gr = genome_gr,
    viral_seqname = "NC_007605.1",
    meta_df = meta_df,
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
