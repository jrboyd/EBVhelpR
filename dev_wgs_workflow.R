library(EBVhelpR)
library(ggplot2)

# Example WGS workflow using EBVhelpR helper functions.
# Update `wgs_root_dir` to your local WGS analysis directory if needed.

wgs_root_dir <- "C:/Users/boydj/OneDrive - UVM Larner College of Medicine/projects_ashley/EBV_DLBCL/P2_viral_WGS"
stopifnot(dir.exists(wgs_root_dir))
# wgs_bwa_dir <- file.path(wgs_root_dir, "output_bwa")
count_dir <- file.path(wgs_root_dir, "chr_read_counts")
bigwig_dir <- file.path(wgs_root_dir, "bigwigs")

# stopifnot(dir.exists(wgs_bwa_dir))
stopifnot(dir.exists(count_dir))
stopifnot(dir.exists(bigwig_dir))

meta_df <- load_meta_data()

# 1) Build reference GRanges from BAM header.
if(FALSE){
    sizes_files = file.path(wgs_root_dir, "chrSizes.txt")
    file.copy(sizes_files, "inst/extdata/WGS_chrSizes.txt")
}

genome_gr <- load_wgs_reference_genome()

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
