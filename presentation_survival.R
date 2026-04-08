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

#version 3 changes
#no combo plots
#phenocycler : just ebna/lmp variants
# viral genome position with highlights for EBER, LMP1, EBNA2, EBNA3

res_dir = "output_survival_v1"
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


meta_df <- load_meta_data()
colors_EBER_status = get_colors_EBER_status()
theme_set(ggpubr::theme_pubr())

#### survival setup ####
# install.packages("survminer")
library(survminer)
library(survival)
surv_df = openxlsx::read.xlsx("inst/extdata/survival_040626.xlsx")

meta_df %>% head
surv_df %>% head
surv_df$Days.until.Death = as.numeric(surv_df$Days.until.Death)
# surv_df = surv_df %>% filter(!is.na(Days.until.Death))

mine_df = surv_df %>% mutate(sample_id = gsub("-", "_", Research.Case.ID)) %>%
    select(sample_id, time = `Days.until.Death`, status = Death.Status)
mine_df = merge(mine_df, meta_df %>% select(sample_id, EBER_status), all.x = TRUE)
mine_df = mine_df %>% filter(EBER_status %in% c("Negative", "Positive"))
fit = survfit(Surv(time, status) ~ EBER_status , data = mine_df)

colors_EBER_status.surv = unlist(colors_EBER_status)
names(colors_EBER_status.surv) = paste0("EBER_status=", names(colors_EBER_status.surv))

p_standard_EBER = ggsurvplot(
    fit,
    data = mine_df,
    size = 1,                 # change line size
    palette = colors_EBER_status.surv, # custom color palettes
    conf.int = TRUE,          # Add confidence interval
    pval = TRUE,               # Add p-value
    # risk.table = TRUE,        # Add risk table
    # risk.table.col = "strata",# Risk table color by groups
    # # legend.labs =
    #     # c("Male", "Female"),    # Change legend labels
    # risk.table.height = 0.25, # Useful to change when you have multiple groups
    ggtheme = theme_bw()      # Change ggplot2 theme
)
p_standard_EBER
false_neg_ids = c(
    "D_EB_18",
    "D_EB_17",
    "D_EB_29",
    "D_EB_21",
    "D_EB_23"
)


mine_df = mine_df %>% mutate(EBER_sensitive = ifelse(sample_id %in% false_neg_ids, "Positive", EBER_status))
# check they flipped


fit.sensitive = survfit(Surv(time, status) ~ EBER_sensitive , data = mine_df)

colors_EBER_status.sens = unlist(colors_EBER_status)
names(colors_EBER_status.sens) = paste0("EBER_sensitive=", names(colors_EBER_status.sens))

p_sensitive_EBER = ggsurvplot(
    fit.sensitive,
    data = mine_df,
    size = 1,                 # change line size
    palette = colors_EBER_status.sens, # custom color palettes
    conf.int = TRUE,          # Add confidence interval
    pval = TRUE,              # Add p-value
    # risk.table = TRUE,        # Add risk table
    # risk.table.col = "strata",# Risk table color by groups
    # # legend.labs =
    #     # c("Male", "Female"),    # Change legend labels
    # risk.table.height = 0.25, # Useful to change when you have multiple groups
    ggtheme = theme_bw()      # Change ggplot2 theme
)
p_sensitive_EBER$plot


p1 = p_standard_EBER$plot + labs(title = "Standard EBER status") +
    scale_fill_manual(labels = function(x)sub(".+=", "", x), values = colors_EBER_status.surv) +
    scale_color_manual(labels = function(x)sub(".+=", "", x), values = colors_EBER_status.surv) +
    labs(x = "Days since Biopsy/Diagnosis")

p2 = p_sensitive_EBER$plot + labs(title = "Sensitive EBER status") +
    scale_fill_manual(labels = function(x)sub(".+=", "", x), values = colors_EBER_status.sens) +
    scale_color_manual(labels = function(x)sub(".+=", "", x), values = colors_EBER_status.sens) +
    labs(x = "Days since Biopsy/Diagnosis")


pg = cowplot::plot_grid(
    nrow = 1,
    p1 + coord_cartesian(xlim = c(0, 4e3)),
    p2 + coord_cartesian(xlim = c(0, 4e3))
)
my_saveplot(pg, "survival", width = 9.5, height = 5)

mine_df %>% filter(sample_id %in% false_neg_ids)

mine_df$sample_id = factor()
mine_df %>% filter(sample_id %in% false_neg_ids) %>% arrange(sample_id)
