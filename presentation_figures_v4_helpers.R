pf_build_subset_dir_name <- function(subset_unique_ids = NULL, run_label = NULL) {
    if (!is.null(run_label) && nzchar(run_label)) {
        clean <- gsub("[^A-Za-z0-9_-]+", "_", run_label)
        return(paste0("run_", clean))
    }

    if (is.null(subset_unique_ids) || !length(subset_unique_ids)) {
        return("all_ids")
    }

    ids <- sort(unique(as.character(subset_unique_ids)))
    preview <- paste(utils::head(ids, 3), collapse = "_")
    preview <- gsub("[^A-Za-z0-9_-]+", "_", preview)
    paste0("subset_n", length(ids), "_", preview)
}

pf_new_context <- function(
        res_dir,
        subset_unique_ids = NULL,
        run_label = NULL,
        wgs_root_dir = NULL,
        cache_dir = "presentation_cache"
) {
    ctx <- new.env(parent = emptyenv())

    ctx$res_dir <- res_dir
    ctx$subset_unique_ids <- if (is.null(subset_unique_ids)) NULL else unique(as.character(subset_unique_ids))
    ctx$run_label <- run_label
    ctx$out_subdir <- pf_build_subset_dir_name(ctx$subset_unique_ids, run_label = ctx$run_label)
    ctx$out_dir <- file.path(res_dir, ctx$out_subdir)
    dir.create(ctx$out_dir, recursive = TRUE, showWarnings = FALSE)

    ctx$wgs_root_dir <- wgs_root_dir
    ctx$cache_dir <- cache_dir
    dir.create(ctx$cache_dir, recursive = TRUE, showWarnings = FALSE)

    ctx$plot_group <- "00"
    ctx$colors_EBER_status <- get_colors_EBER_status()
    ctx$meta_df <- load_meta_data() %>% filter(sample_id %in% subset_unique_ids)

    ggplot2::theme_set(ggpubr::theme_pubr())

    ctx
}

pf_res_file <- function(ctx, f) {
    f <- file.path(dirname(f), paste0(ctx$plot_group, "_", f))
    out_f <- file.path(ctx$out_dir, f)
    dir.create(dirname(out_f), recursive = TRUE, showWarnings = FALSE)
    out_f
}

pf_increase_plot_group <- function(ctx) {
    next_group <- format(as.numeric(ctx$plot_group) + 1, digits = 2, width = 2)
    ctx$plot_group <- gsub(" ", "0", next_group)
    invisible(ctx$plot_group)
}

pf_saveplot <- function(ctx, plot, name, width, height) {
    ggplot2::ggsave(pf_res_file(ctx, paste0(name, ".png")), plot, width = width, height = height)
    ggplot2::ggsave(pf_res_file(ctx, paste0(name, ".pdf")), plot, width = width, height = height)
}

pf_select_scope_ids <- function(cq, valid_statuses, valid_types) {
    sel <- cq@summary_df %>%
        subset(probe_control == "" & EBER_status %in% valid_statuses & sample_type %in% valid_types)
    unique(as.character(sel$unique_id))
}

pf_section_wgs_setup <- function(ctx) {
    ctx$wgs_files_df <- setup_wgs_files(
        wgs_root_dir = ctx$wgs_root_dir,
        meta_df = ctx$meta_df
    )

    #filter out not performed EBER_status
    ctx$wgs_files_df <- ctx$wgs_files_df %>%
        dplyr::filter(EBER_status %in% c("Positive", "Negative"))

    ctx$genome_gr <- load_wgs_reference_genome()

    viral_genes_gr <- rtracklayer::import.gff("inst/extdata/NC_007605.1_gene_reference.gff3")
    viral_genes_df <- viral_genes_gr %>%
        as.data.frame() %>%
        dplyr::mutate(gene = ifelse(grepl("EBER", product), sub(" .+", "", product), gene))

    tmp1 <- viral_genes_df %>% dplyr::filter(grepl("EBER", product) & type == "transcript")
    tmp2 <- viral_genes_df %>% dplyr::filter((grepl("EBNA.+[0-9].?", gene) | grepl("LMP", gene)) & type == "exon")
    viral_genes_df.highlight <- rbind(tmp1, tmp2) %>%
        dplyr::mutate(gene = ifelse(grepl("EBER", product), sub(" .+", "", product), gene))

    viral_genes_label <- viral_genes_df.highlight %>%
        dplyr::group_by(gene) %>%
        dplyr::summarise(x = mean(start + end) / 2, .groups = "drop")

    ctx$viral_genes_df = viral_genes_df
    ctx$viral_genes_df.highlight = viral_genes_df.highlight
    ctx$viral_genes_plot = plot_viral_gene_ref(ctx$viral_genes_df, ctx$viral_genes_df.highlight$gene)

    pf_increase_plot_group(ctx)
    pf_saveplot(ctx, ctx$viral_genes_plot, "ebv_reference_highlights", width = 8.15, height = 2.65)

    ctx$wgs_count_summary <- load_wgs_count_summary(
        wgs_files_df = ctx$wgs_files_df,
        genome_gr = ctx$genome_gr,
        viral_seqname = "NC_007605.1"
    )

    pileup_cache <- file.path(ctx$cache_dir, "pileup_df.Rds")
    if (file.exists(pileup_cache)) {
        ctx$pileup_df <- readRDS(pileup_cache)
    } else {
        ctx$pileup_df <- load_wgs_bigwig_pileup(
            wgs_files_df = ctx$wgs_files_df,
            genome_gr = ctx$genome_gr,
            viral_seqname = "NC_007605.1",
            smooth_n = 50,
            mc_cores = 4
        )
        saveRDS(ctx$pileup_df, pileup_cache)
    }

    stopifnot(is.factor(ctx$wgs_count_summary$sample_id))
    ctx$id_levels = levels(ctx$wgs_count_summary$sample_id)
    invisible(ctx)
}

append_viral_reference = function(in_plot,
                                  ctx,
                                  rel_heights,
                                  xrng,
                                  yrng){
    plist = seqsetvis:::sync_height(list(
        in_plot + coord_cartesian(xlim = xrng, ylim = yrng) +
            ggplot2::scale_x_continuous(labels = function(x) x / 1000) +
            ggplot2::labs(x = "EBV genomic position (kb)") +
            labs(x = "") + theme(axis.text.x = element_blank()),
        ctx$viral_genes_plot + coord_cartesian(xlim = xrng)
    ), sync_width = TRUE)
    cowplot::plot_grid(plotlist = plist, ncol = 1, rel_heights = rel_heights)

}

.plot_wgs = function(ctx, wgs_summary, wgs_pileup, prefix){
    p_enrichment <- plot_wgs_viral_enrichment(wgs_summary) +
        ggplot2::scale_fill_manual(values = ctx$colors_EBER_status) +
        ggplot2::scale_color_manual(values = ctx$colors_EBER_status) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(hjust = 0, size = 8))

    p_lines <- plot_wgs_pileup_lines(wgs_pileup, ylim = c(0, 2000))

    xrng = range(wgs_pileup$x)

    pg_lines.2k = append_viral_reference(p_lines, ctx, c(2, 1.5), xrng = xrng, yrng = c(0, 2000))
    pg_lines.zoom = append_viral_reference(p_lines, ctx, c(2, 1.5), xrng = xrng, yrng = c(0, 20))



    pg_lines <- cowplot::plot_grid(pg_lines.2k, pg_lines.zoom, nrow = 1)

    p_heat <- plot_wgs_pileup_heatmap(
        pileup_df = wgs_pileup,
        wgs_count_summary = wgs_summary,
        normalize_by_sample = FALSE,
        status_colors = ctx$colors_EBER_status
    ) +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::theme(axis.text.y = ggplot2::element_text(hjust = 0, size = 8))

    p_heat_norm <- plot_wgs_pileup_heatmap(
        pileup_df = wgs_pileup,
        wgs_count_summary = wgs_summary,
        normalize_by_sample = TRUE,
        status_colors = ctx$colors_EBER_status
    ) +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::theme(axis.text.y = ggplot2::element_text(hjust = 0, size = 8))

    p_heat = append_viral_reference(p_heat, ctx, rel_heights = c(3, 1), xrng = xrng, yrng = NULL)
    p_heat_norm = append_viral_reference(p_heat_norm, ctx, rel_heights = c(3, 1), xrng = xrng, yrng = NULL)

    pf_increase_plot_group(ctx)
    pf_saveplot(ctx, p_enrichment, paste0(prefix, "_wgs_enrichment_barplot"), width = 11, height = 8)
    pf_saveplot(ctx, pg_lines, paste0(prefix, "_wgs_profiles_lineplot"), width = 12.7, height = 6.2)
    pf_saveplot(ctx, p_heat, paste0(prefix, "_wgs_profiles_heatmap"), width = 11, height = 8+3)
    pf_saveplot(ctx, p_heat_norm, paste0(prefix, "_wgs_profiles_heatmap_normalized"), width = 11, height = 8+3)

    invisible(ctx)
}

pf_section_wgs_initial_plots <- function(ctx) {
    invisible(.plot_wgs(ctx, ctx$wgs_count_summary, ctx$pileup_df, "initial"))
}

pf_section_overlap_selection <- function(ctx) {
    ctx$cq_r4 <- CellQuery(assay_type = EBV_ASSAY_TYPES$RNAScope_4plex)
    ids_r4 <- pf_select_scope_ids(ctx$cq_r4, ctx$meta_df$EBER_status, ctx$meta_df$sample_type)
    ctx$cq_r4 <- set_selected_unique_ids(ctx$cq_r4, ctx$id_levels)


    ctx$cq_r3i <- CellQuery(assay_type = EBV_ASSAY_TYPES$`RNAScope_3plex+IF`)
    ids_r3i <- pf_select_scope_ids(ctx$cq_r3i, ctx$wgs_count_summary$EBER_status, ctx$meta_df$sample_type)
    ctx$cq_r3i <- set_selected_unique_ids(ctx$cq_r3i, ctx$id_levels)

    ctx$cq_p <- CellQuery(assay_type = EBV_ASSAY_TYPES$Phenocycler)
    ids_p <- pf_select_scope_ids(ctx$cq_p, ctx$wgs_count_summary$EBER_status, ctx$meta_df$sample_type)
    ctx$cq_p <- set_selected_unique_ids(ctx$cq_p, ctx$id_levels)

    all_ids <- list(
        WGS = ctx$id_levels,
        `RNAscope 4plex` = ids_r4,
        `RNAscope 3plexIF` = ids_r3i,
        Phenocycler = ids_p
    )

    pf_increase_plot_group(ctx)
    upr <- seqsetvis::ssvFeatureUpset(all_ids, return_UpSetR = TRUE)

    grDevices::pdf(pf_res_file(ctx, "assay_overlap_upset.pdf"), width = 6, height = 5)
    print(upr)
    grDevices::dev.off()

    grDevices::png(pf_res_file(ctx, "assay_overlap_upset.png"), width = 6, height = 5, units = "in", res = 200)
    print(upr)
    grDevices::dev.off()

    membs <- seqsetvis::ssvMakeMembTable(all_ids)
    in_all <- apply(membs, 1, all)
    final_ids <- rownames(membs)[in_all]

    if (!is.null(ctx$subset_unique_ids) && length(ctx$subset_unique_ids)) {
        requested <- unique(as.character(ctx$subset_unique_ids))
        missing_ids <- setdiff(requested, final_ids)
        if (length(missing_ids)) {
            warning(
                "Requested subset unique_id values not found in overlap set: ",
                paste(missing_ids, collapse = ", ")
            )
        }
        final_ids <- intersect(final_ids, requested)
    }

    ctx$final_ids <- final_ids
    saveRDS(ctx$final_ids, pf_res_file(ctx, "final_ids.Rds"))

    invisible(ctx)
}

pf_section_apply_final_ids <- function(ctx) {
    ctx$wgs_count_summary.re <- ctx$wgs_count_summary %>%
        dplyr::filter(sample_id %in% ctx$final_ids)
    ctx$wgs_count_summary.re$sample_id <- droplevels(ctx$wgs_count_summary.re$sample_id)

    ctx$pileup_df.re <- ctx$pileup_df %>%
        dplyr::filter(sample %in% ctx$final_ids)

    ctx$cq_r4 <- set_selected_unique_ids(ctx$cq_r4, ctx$final_ids)
    ctx$cq_r3i <- set_selected_unique_ids(ctx$cq_r3i, ctx$final_ids)
    ctx$cq_p <- set_selected_unique_ids(ctx$cq_p, ctx$final_ids)

    stopifnot(is.factor(ctx$wgs_count_summary.re$sample_id))
    ctx$id_levels = levels(ctx$wgs_count_summary.re$sample_id)

    invisible(ctx)
}

pf_section_wgs_replots <- function(ctx) {
    invisible(.plot_wgs(ctx, ctx$wgs_count_summary.re, ctx$pileup_df.re, "restricted"))
}

pf_plot_scope_summary_bars <- function(cq, id_lev, status_colors) {
    qsum <- get_query_summary_df(cq)
    qsum$sample_id <- factor(qsum$sample_id, levels = id_lev)
    qsum <- qsum %>% subset(!grepl("_", combo))

    ggplot2::ggplot(qsum, ggplot2::aes(x = Positive_Percent, y = sample_id, fill = EBER_status)) +
        ggplot2::geom_col() +
        ggplot2::facet_wrap(~combo, scales = "free_x", nrow = 1) +
        ggplot2::scale_fill_manual(values = status_colors) +
        ggplot2::labs(x = "Percent Positive", y = "", fill = "EBER status")
}

pf_plot_scope_summary_boxplot <- function(cq, status_colors) {
    qsum <- get_query_summary_df(cq)
    qsum <- qsum %>% subset(!grepl("_", combo))

    ggplot2::ggplot(qsum, ggplot2::aes(x = combo, y = Positive_Percent, fill = EBER_status)) +
        ggplot2::geom_boxplot() +
        ggplot2::scale_fill_manual(values = status_colors) +
        ggplot2::labs(x = "Probe", y = "Percent Positive", fill = "EBER status")
}

pf_section_scope_summary_plots <- function(ctx) {
    p_sum_bars_r4 <- pf_plot_scope_summary_bars(ctx$cq_r4, ctx$id_levels, ctx$colors_EBER_status)
    p_sum_bars_r3i <- pf_plot_scope_summary_bars(ctx$cq_r3i, ctx$id_levels, ctx$colors_EBER_status)

    p_sum_box_r4 <- pf_plot_scope_summary_boxplot(ctx$cq_r4, ctx$colors_EBER_status)
    p_sum_box_r3i <- pf_plot_scope_summary_boxplot(ctx$cq_r3i, ctx$colors_EBER_status)

    qsum_p <- get_query_summary_df(ctx$cq_p)
    pos_p <- qsum_p %>%
        dplyr::select(sample_id, dplyr::contains("Cells") & dplyr::contains("%")) %>%
        tidyr::pivot_longer(cols = !sample_id) %>%
        dplyr::mutate(name = sub("% ", "", name)) %>%
        dplyr::mutate(name = sub(" Positive Cells", "", name))

    pos_p <- merge(ctx$meta_df, pos_p)
    pos_p$sample_id <- factor(pos_p$sample_id, levels = ctx$id_levels)

    pos_p.ebv = subset(pos_p %>% filter(grepl("(LMP)|(EBN)", name)))

    p_sum_bars_p <- ggplot2::ggplot(pos_p, ggplot2::aes(x = value, y = sample_id, fill = EBER_status)) +
        ggplot2::geom_col() +
        ggplot2::facet_wrap(~name, scales = "free_x", nrow = 2) +
        ggplot2::scale_fill_manual(values = ctx$colors_EBER_status) +
        ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(2)) +
        ggplot2::labs(x = "Percent Positive", y = "", fill = "EBER status") +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8), axis.text.x = ggplot2::element_text(size = 8))

    p_sum_bars_p.ebv = ggplot2::ggplot(pos_p.ebv, ggplot2::aes(x = value, y = sample_id, fill = EBER_status)) +
        ggplot2::geom_col() +
        ggplot2::facet_wrap(~name, scales = "free_x", nrow = 2) +
        ggplot2::scale_fill_manual(values = ctx$colors_EBER_status) +
        ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(2)) +
        ggplot2::labs(x = "Percent Positive", y = "", fill = "EBER status") +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8), axis.text.x = ggplot2::element_text(size = 8))

    p_sum_box_p <- ggplot2::ggplot(pos_p, ggplot2::aes(x = name, y = value, fill = EBER_status)) +
        ggplot2::geom_boxplot() +
        ggplot2::scale_fill_manual(values = ctx$colors_EBER_status) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)) +
        ggplot2::labs(y = "Percent Positive", x = "Probe", fill = "EBER status")

    p_sum_box_p.ebv = ggplot2::ggplot(pos_p.ebv, ggplot2::aes(x = name, y = value, fill = EBER_status)) +
        ggplot2::geom_boxplot() +
        ggplot2::scale_fill_manual(values = ctx$colors_EBER_status) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)) +
        ggplot2::labs(y = "Percent Positive", x = "Probe", fill = "EBER status")

    pf_increase_plot_group(ctx)
    pf_saveplot(ctx, p_sum_bars_r4, name = "RNAScope_4plex_barplot", width = 11, height = 7.2)
    pf_saveplot(ctx, p_sum_bars_r3i, name = "RNAScope_3plexIF_barplot", width = 11, height = 7.2)
    pf_saveplot(ctx, p_sum_bars_p, name = "Phenocycler_barplot.all_probes", width = 11, height = 8.3)
    pf_saveplot(ctx, p_sum_bars_p.ebv, name = "Phenocycler_barplot", width = 5.4, height = 8.3)

    pf_saveplot(ctx, p_sum_box_r4, name = "RNAScope_4plex_boxplot", width = 8.4, height = 5)
    pf_saveplot(ctx, p_sum_box_r3i, name = "RNAScope_3plexIF_boxplot", width = 8.4, height = 5)
    pf_saveplot(ctx, p_sum_box_p, name = "Phenocycler_boxplot.all_probes", width = 7.5, height = 6.4)
    pf_saveplot(ctx, p_sum_box_p.ebv, name = "Phenocycler_boxplot", width = 4, height = 4.5)
    pf_saveplot(ctx, p_sum_box_p.ebv + coord_cartesian(ylim = c(0, 15)), name = "Phenocycler_boxplot.zoom", width = 4, height = 4.5)

    invisible(ctx)
}

pf_plot_enrichment_vs_wgs <- function(df, grouped = FALSE, ncol = NULL, nrow = NULL, drop = TRUE) {
    aes_obj <- if (grouped) {
        ggplot2::aes(x = viral_enrichment, y = Positive_Percent, color = EBER_status, group = combo)
    } else {
        ggplot2::aes(x = viral_enrichment, y = Positive_Percent, color = EBER_status)
    }

    ggplot2::ggplot(df, aes_obj) +
        ggplot2::geom_point() +
        ggplot2::facet_wrap(~combo, scales = "free_y", ncol = ncol, nrow = nrow, drop = drop) +
        ggplot2::scale_x_log10() +
        ggplot2::scale_y_log10() +
        ggpmisc::stat_poly_line(method = "lm", formula = y ~ x) +
        ggpmisc::stat_poly_eq(
            ggplot2::aes(label = ..rr.label..),
            formula = y ~ x,
            parse = TRUE,
            label.y = c(.8, .95)
        )
}

pf_section_save_cell_queries <- function(ctx) {
    pf_increase_plot_group(ctx)
    saveRDS(ctx$cq_r4, pf_res_file(ctx, "cq_r4.Rds"))
    saveRDS(ctx$cq_r3i, pf_res_file(ctx, "cq_r3i.Rds"))
    saveRDS(ctx$cq_p, pf_res_file(ctx, "cq_p.Rds"))
    invisible(ctx)
}

pf_section_scope_correlations <- function(ctx) {
    pf_increase_plot_group(ctx)

    qsum_r4 <- get_query_summary_df(ctx$cq_r4)
    xy_df <- merge(
        qsum_r4,
        ctx$wgs_count_summary.re %>% dplyr::select(sample_id, viral_enrichment),
        by = "sample_id"
    )
    xy_df$sample_id <- factor(xy_df$sample_id, levels = ctx$id_levels)
    xy_df <- xy_df %>% dplyr::filter(!grepl("Unstai", combo))

    p_grouped <- pf_plot_enrichment_vs_wgs(xy_df, grouped = FALSE)
    p_ungrouped <- pf_plot_enrichment_vs_wgs(xy_df, grouped = TRUE)
    pf_saveplot(ctx, p_grouped, name = "RNAscope4plex_vs_WGS_grouped", width = 9.3, height = 7.5)
    pf_saveplot(ctx, p_ungrouped, name = "RNAscope4plex_vs_WGS_ungrouped", width = 9.3, height = 7.5)

    qsum_r3i <- get_query_summary_df(ctx$cq_r3i)
    xy_df <- merge(
        qsum_r3i,
        ctx$wgs_count_summary.re %>% dplyr::select(sample_id, viral_enrichment),
        by = "sample_id"
    )
    xy_df$sample_id <- factor(xy_df$sample_id, levels = ctx$id_levels)
    xy_df <- xy_df %>% dplyr::filter(!grepl("Unstai", combo))

    p_grouped <- pf_plot_enrichment_vs_wgs(xy_df, grouped = FALSE)
    p_ungrouped <- pf_plot_enrichment_vs_wgs(xy_df, grouped = TRUE)
    pf_saveplot(ctx, p_grouped, name = "RNAscope3plexIF_vs_WGS_grouped", width = 9.3, height = 7.5)
    pf_saveplot(ctx, p_ungrouped, name = "RNAscope3plexIF_vs_WGS_ungrouped", width = 9.3, height = 7.5)

    qsum_p <- get_query_summary_df(ctx$cq_p)
    pos_p <- qsum_p %>%
        dplyr::select(sample_id, dplyr::contains("Cells") & dplyr::contains("%")) %>%
        tidyr::pivot_longer(cols = !sample_id) %>%
        dplyr::mutate(name = sub("% ", "", name)) %>%
        dplyr::mutate(name = sub(" Positive Cells", "", name))

    pos_p <- merge(ctx$meta_df, pos_p)
    pos_p$sample_id <- factor(pos_p$sample_id, levels = ctx$id_levels)

    qsum <- pos_p %>% dplyr::rename(Positive_Percent = value, combo = name)
    xy_df <- merge(
        qsum,
        ctx$wgs_count_summary.re %>% dplyr::select(sample_id, viral_enrichment),
        by = "sample_id"
    )
    xy_df$sample_id <- factor(xy_df$sample_id, levels = ctx$id_levels)
    xy_df <- xy_df %>% dplyr::filter(!grepl("Unstai", combo))

    pheno_combos <- unique(xy_df$combo)
    n_per <- 9
    todo <- split(pheno_combos, ceiling(seq_along(pheno_combos) / n_per))
    xy_df_full <- xy_df

    for (i in seq_along(todo)) {
        sel_combos <- todo[[i]]
        xy_df <- xy_df_full %>% dplyr::filter(combo %in% sel_combos)
        n_combos <- length(unique(xy_df$combo))

        if (n_combos < 9) {
            n_add <- 9 - n_combos
            xy_df$combo <- factor(xy_df$combo, levels = c(unique(xy_df$combo), paste("empty", LETTERS[seq(n_add)])))
        }

        f_suffix <- paste0(".", i, "_of_", length(todo))
        p_grouped <- pf_plot_enrichment_vs_wgs(xy_df, grouped = FALSE, ncol = 3, nrow = 3, drop = FALSE)
        p_ungrouped <- pf_plot_enrichment_vs_wgs(xy_df, grouped = TRUE, ncol = 3, nrow = 3, drop = FALSE)
        pf_saveplot(ctx, p_grouped, name = paste0("Phenocycler_vs_WGS_grouped", f_suffix), width = 9.3, height = 7.5)
        pf_saveplot(ctx, p_ungrouped, name = paste0("Phenocycler_vs_WGS_ungrouped", f_suffix), width = 9.3, height = 7.5)
    }

    invisible(ctx)
}

pf_run_all_sections <- function(ctx) {
    ctx = pf_section_wgs_setup(ctx)
    ctx = pf_section_wgs_initial_plots(ctx)
    ctx = pf_section_overlap_selection(ctx)
    ctx = pf_section_apply_final_ids(ctx)
    ctx = pf_section_wgs_replots(ctx)
    ctx = pf_section_scope_summary_plots(ctx)
    ctx = pf_section_save_cell_queries(ctx)
    ctx = pf_section_scope_correlations(ctx)
    invisible(ctx)
}

pf_run_sections <- function(ctx, sections) {
    runners <- list(
        wgs_setup = pf_section_wgs_setup,
        wgs_initial_plots = pf_section_wgs_initial_plots,
        overlap_selection = pf_section_overlap_selection,
        apply_final_ids = pf_section_apply_final_ids,
        wgs_replots = pf_section_wgs_replots,
        scope_summary_plots = pf_section_scope_summary_plots,
        save_cell_queries = pf_section_save_cell_queries,
        scope_correlations = pf_section_scope_correlations
    )

    for (nm in sections) {
        if (!nm %in% names(runners)) {
            stop("Unknown section: ", nm)
        }
        ctx = runners[[nm]](ctx)
    }

    invisible(ctx)
}

plot_viral_gene_ref = function(viral_genes_df, highlight_genes){
    viral_genes_df = viral_genes_df %>% filter(type == "exon")
    viral_genes_df = viral_genes_df %>% mutate(gene_label = ifelse(gene %in% highlight_genes, gene, "other"))
    head(viral_genes_df)
    viral_genes_df$gene_label %>% table
    viral_genes_df$gene_label = factor(viral_genes_df$gene_label)
    ?relevel

    lev_o =c(setdiff(levels(viral_genes_df$gene_label), "other"), "other")
    lev_o = rev(lev_o)
    viral_genes_df$gene_label = factor(viral_genes_df$gene_label, levels = lev_o)
    viral_genes_df$gene_label
    ggplot(viral_genes_df, aes(x = start, xend = end, y = gene_label, yend = gene_label)) +
        geom_segment(linewidth = 2)

    ggplot(viral_genes_df, aes(x = start, xend = end, y = gene_label, yend = gene_label)) +
        geom_point()

    rect_height = .8
    viral_genes_df = viral_genes_df %>% mutate(ymin = as.numeric(gene_label)-rect_height/2, ymax = as.numeric(gene_label)+rect_height/2)

    gene_names = levels(viral_genes_df$gene_label)
    gene_names = split(gene_names, sub("-.+", "", gene_names))
    gene_colors = gene_names
    gene_colors$EBER = seqsetvis::safeBrew(c(1:5, gene_names$EBER), pal = "Blues")
    gene_colors$EBNA = seqsetvis::safeBrew(c(1:3, gene_names$EBNA), pal = "Reds")
    gene_colors$LMP = seqsetvis::safeBrew(c(1:4, gene_names$LMP), pal = "Greens")
    gene_colors$other = c("other" = "gray")
    names(gene_colors) = NULL
    gene_colors = unlist(gene_colors)

    ggplot(viral_genes_df, aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax)) +
        geom_point(aes(x = start, y = gene_label), alpha = 0) +
        geom_rect(linewidth = .5, aes(color = gene_label, fill = gene_label), show.legend = FALSE) +
        scale_color_manual(values = gene_colors) +
        scale_fill_manual(values = gene_colors) +
        labs(y = "", x = "") +
        ggplot2::scale_x_continuous(labels = function(x) x / 1000) +
        ggplot2::labs(x = "EBV genomic position (kb)")

}
