library(EBVhelpR)
library(ggplot2)
library(tidyverse)

source("presentation_figures_v4_helpers.R")

# Change this to your preferred output root.
res_dir <- "output_presentation_04092026_v6"

meta_df = load_meta_data()
meta_df = meta_df %>% filter(sample_type != "control")

all_ids = meta_df$sample_id
ids_by_type = split(all_ids, sub("_.+", "", all_ids))
run_names = c("CTEBV" = "control", "D" = "DLBCL")
names(ids_by_type) = run_names[names(ids_by_type)]

ids_by_type$all_samples = all_ids

name = names(ids_by_type)[2]
name = "control"
# run controls and DLBCL separate
for(name in names(ids_by_type)){
    message("running ", name)
    ctx <- pf_new_context(
        res_dir = res_dir,
        subset_unique_ids = ids_by_type[[name]],
        run_label = name
    )

    ctx = pf_section_wgs_setup(ctx)
    ctx = pf_section_wgs_initial_plots(ctx)
    ctx = pf_section_overlap_selection(ctx)
    ctx = pf_section_apply_final_ids(ctx)
    ctx = pf_section_wgs_replots(ctx)
    ctx = pf_section_scope_summary_plots(ctx)
    ctx = pf_section_save_cell_queries(ctx)
    ctx = pf_section_scope_correlations(ctx)

    # Run complete workflow.
    # pf_run_all_sections(ctx)
    message("Outputs written to: ", ctx$out_dir)
}



# Or run selected sections, in order, for partial reruns.
# pf_run_sections(ctx, c(
#     "wgs_setup",
#     "wgs_initial_plots",
#     "overlap_selection",
#     "apply_final_ids",
#     "wgs_replots",
#     "scope_summary_plots",
#     "save_cell_queries",
#     "scope_correlations"
# ))

ctx$pileup_df

viral_genes_label <- viral_genes_df.highlight %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(x = mean(start + end) / 2, .groups = "drop")

p_ref <- ggplot2::ggplot(viral_genes_df) +
    ggplot2::geom_segment(ggplot2::aes(x = start, xend = end, y = 0, yend = 0), linewidth = 1.2) +
    ggplot2::geom_segment(data = viral_genes_df.highlight, ggplot2::aes(x = start, xend = end, y = 0, yend = 0), color = "red", linewidth = 3) +
    ggrepel::geom_label_repel(data = viral_genes_label, ggplot2::aes(x = x, y = 0, label = gene))
p_ref


high_df = ctx$viral_genes_df
high_df = high_df %>% filter(type == "exon")
high_df = high_df %>% mutate(gene_label = ifelse(gene %in% ctx$viral_genes_df.highlight$gene, gene, "other"))
head(high_df)
high_df$gene_label %>% table
high_df$gene_label = factor(high_df$gene_label)
?relevel

lev_o =c(setdiff(levels(high_df$gene_label), "other"), "other")
lev_o = rev(lev_o)
high_df$gene_label = factor(high_df$gene_label, levels = lev_o)
high_df$gene_label
ggplot(high_df, aes(x = start, xend = end, y = gene_label, yend = gene_label)) +
    geom_segment(linewidth = 2)

ggplot(high_df, aes(x = start, xend = end, y = gene_label, yend = gene_label)) +
    geom_point()

rect_height = .8
high_df = high_df %>% mutate(ymin = as.numeric(gene_label)-rect_height/2, ymax = as.numeric(gene_label)+rect_height/2)

gene_names = levels(high_df$gene_label)
gene_names = split(gene_names, sub("-.+", "", gene_names))
gene_colors = gene_names
gene_colors$EBER = seqsetvis::safeBrew(c(1:5, gene_names$EBER), pal = "Blues")
gene_colors$EBNA = seqsetvis::safeBrew(c(1:3, gene_names$EBNA), pal = "Reds")
gene_colors$LMP = seqsetvis::safeBrew(c(1:4, gene_names$LMP), pal = "Greens")
gene_colors$other = c("other" = "gray")
names(gene_colors) = NULL
gene_colors = unlist(gene_colors)

ggplot(high_df, aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax)) +
    geom_point(aes(x = start, y = gene_label), alpha = 0) +
    geom_rect(linewidth = .5, aes(color = gene_label, fill = gene_label), show.legend = FALSE) +
    scale_color_manual(values = gene_colors) +
    scale_fill_manual(values = gene_colors) +
    labs(y = "", x = "") +
    ggplot2::scale_x_continuous(labels = function(x) x / 1000) +
    ggplot2::labs(x = "EBV genomic position (kb)")

####


