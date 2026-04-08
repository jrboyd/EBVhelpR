library(EBVhelpR)
library(ggplot2)
library(tidyverse)

source("presentation_figures_v4_helpers.R")

# Change this to your preferred output root.
res_dir <- "output_presentation_04092026_v4"

meta_df = load_meta_data()
meta_df = meta_df %>% filter(sample_type != "control")

all_ids = meta_df$sample_id
ids_by_type = split(all_ids, sub("_.+", "", all_ids))
run_names = c("CTEBV" = "control", "D" = "DLBCL")
names(ids_by_type) = run_names[names(ids_by_type)]

ids_by_type$all_samples = all_ids

name = names(ids_by_type)[2]
# run controls and DLBCL separate
for(name in names(ids_by_type)){
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


