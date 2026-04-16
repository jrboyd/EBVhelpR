library(EBVhelpR)
library(ggplot2)
library(tidyverse)

res_dir = "output_survival_v2"
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
    ggsave(res_file(paste0(name, ".png")), plot, width = width, height = height, bg = "white")
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
dim(surv_df)

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

mine_df %>% filter(sample_id %in% false_neg_ids) %>% arrange(sample_id)

days_to_test = c(500, 1000, 2000)
alive_at_days = 500

eber_to_test = c("EBER_status", "EBER_sensitive")
eber_var = "EBER_status"
eber_var = "EBER_sensitive"


tables_list = list()
pvalues_list = list()

for(alive_at_days in days_to_test){
    for(eber_var in eber_to_test){
        #remove unknown status samples
        fish_df = mine_df %>% filter(!is.na(time))
        fish_df = fish_df %>% mutate(incomplete = (time < alive_at_days & status == 0))
        fish_df = fish_df %>% mutate(still_alive = time > alive_at_days)
        cont_df = fish_df %>% group_by(!!sym(eber_var), still_alive) %>% summarise(N = length(sample_id))
        cont_df = cont_df %>% pivot_wider(id_cols = !!sym(eber_var), names_from = still_alive, values_from = N, values_fill = 0)
        cont_mat = as.matrix(cont_df[,-1])
        rownames(cont_mat) = cont_df[[eber_var]]

        table_df = reshape2::melt(cont_mat)
        colnames(table_df) = c("EBER", "alive", "count")

        name = paste(eber_var, alive_at_days)
        tables_list[[name]] = table_df


        test_res = fisher.test(cont_mat, alternative = "two.sided")
        test_res$p.value
        pvalues_list[[name]] = test_res$p.value

    }
}

tables_list
pvalues_list
mat_df = data.frame(res_names = names(tables_list)) %>%
    separate(col = res_names,
             into = c("EBER_type", "days"),
             sep =  " ", remove = FALSE)

table_df = tables_list[[1]]
tab_plots = lapply(names(tables_list), function(name){
    table_df = tables_list[[name]]
    pval = pvalues_list[[name]]
    table_df$EBER = factor(table_df$EBER, levels = rev(levels(table_df$EBER)))

    table_df$alive_label = ifelse(table_df$alive, "Alive", "Deceased")
    table_df$alive_label = factor(table_df$alive_label, levels = c("Deceased", "Alive"))

    ggplot(table_df) +
        geom_text(aes(x = EBER, y = alive_label, label = count)) +
        scale_x_discrete(position = "top") +
        theme(axis.line = element_blank(), axis.ticks = element_blank()) +
        theme(axis.text = element_text(size = 8)) +
        labs(x = "", y = "", caption = paste("p-value", format(pval, digits = 5)))
})

pg_main = cowplot::plot_grid(plotlist = tab_plots, byrow = TRUE, ncol = 2)

row_labels = lapply(days_to_test, function(x){
    ggplot() + annotate("text", x= .7, y = .5, label = x) + theme_void() +
        coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
})

status_rename = c(
    "EBER_status" = "Standard",
    "EBER_sensitive" = "Sensitive"
)
col_labels = lapply(c("", status_rename[eber_to_test]), function(x){
    ggplot() + annotate("text", x= .5, y = .25, label = x) + theme_void() +
        coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
})

pg_row_labels = cowplot::plot_grid(plotlist = row_labels, ncol = 1)
pg_main_row  = cowplot::plot_grid(pg_row_labels, pg_main, nrow = 1, rel_widths = c(1, 5))
pg_col_labels = cowplot::plot_grid(plotlist = col_labels, nrow = 1, rel_widths = c(1, 5))

pg_assembly = cowplot::plot_grid(
    pg_col_labels,
    pg_main_row,
    ncol = 1, rel_heights = c(1, 6)
)

pg_assembly.final = pg_assembly +
    cowplot::draw_text(text = "Days Since Diagnosis/Biopsy", x = .05, y = .42, angle = 90, vjust = 0) +
    cowplot::draw_text(text = "EBER Status Type", x = .6, y = .98, angle = 0, vjust = 1)

my_saveplot(cowplot::as_grob(pg_assembly.final), "fisher_test", width = 5, height = 6.4)
