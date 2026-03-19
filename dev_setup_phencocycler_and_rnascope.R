library(magrittr)
library(tidyverse)

EBVhelpR::write_all_package_data()

EBVhelpR::load_phenocycler_summary_files()
# debug(EBVhelpR::harmonize_phenocycler_summary_files)
EBVhelpR::harmonize_phenocycler_summary_files()

rscop_summary = EBVhelpR::load_rnascope_summary_files()
rscop_summary$source %>% table
EBVhelpR::harmonize_rnascope_summary_files()


debug(load_cell_source_files)
load_cell_source_files()

#' 
#' #' Title
#' #'
#' #' @returns
#' #' @export
#' #'
#' #' @examples
#' load_cell_source_files = function(){
#'   
#'   pkg_data_dir = get_pkg_data_dir()
#'   
#'   all_cell_data_files = list.files(pkg_data_dir, pattern = ".cell_data.csv$", recursive = TRUE)
#'   
#'   cell_meta_files = list.files(pkg_data_dir, pattern = "meta.+csv$", recursive = TRUE, full.names = TRUE)
#'   lapply(cell_meta_files, read.csv)
#'   
#'   
#'   cell_df = data.frame(file = all_cell_data_files)
#'   
#'   cell_df = cell_df %>% separate(file, sep = "/", into = c("assay", "image_name"), remove = FALSE) %>% mutate(image_name = sub("\\..+", "", image_name))
#'   cell_df$name = cell_df$image_name
#'   cell_df = cell_df %>% 
#'     mutate(name = sub("_ ?[cC]ro.+", "", name)) %>% 
#'     mutate(name = sub("^[0-9]+_", "", name)) %>% 
#'     mutate(name = sub("EBER-LMP1-EBNA1", "PosCTL", name)) %>% 
#'     mutate(name = sub("_Rescanned", "", name)) %>% 
#'     mutate(name = gsub("_Rescanned", "", name))
#'   cell_df = cell_df %>% mutate(sample_type = "cohort") %>%
#'     # mutate(sample_type = ifelse(grepl("CellPellet", name), "pellet", sample_type)) %>%
#'     mutate(sample_type = ifelse(grepl("CTL", name), "control", sample_type)) 
#'   
#'   cell_df = cell_df %>% 
#'     mutate(name = gsub("-", "_", name)) %>% 
#'     mutate(name = sub("^CTEBV", "CTEBV_", name))  %>% 
#'     mutate(name = sub("_Scan2", "", name)) %>%
#'     mutate(name = sub("CellPellet_", "", name))
#'   
#'   cell_df$sample_id = cell_df$name
#'   
#'   # intersect(meta_df$sample_id, cell_df$sample_id)
#'   # setdiff(meta_df$sample_id, cell_df$sample_id)
#'   # setdiff(cell_df$sample_id, meta_df$sample_id)
#'   
#'   
#'   split(cell_df, cell_df$assay)
#'   cell_df.by_type = split(cell_df, cell_df$sample_type)
#'   
#'   meta_df = EBVhelpR::load_meta_data()
#'   
#'   cell_df.by_type$cohort = merge(cell_df.by_type$cohort, meta_df, all.x = TRUE)
#'   
#'   #cohort sample are all set
#'   stopifnot(all(cell_df.by_type$cohort$sample_id %in% meta_df$sample_id))
#'   
#'   
#'   
#'   #control samples need additional tweaks
#'   cell_df.by_type$control = cell_df.by_type$control %>% 
#'     mutate(sample_id  = sub("_2$", "", sample_id)) %>% 
#'     mutate(control_type = ifelse(grepl("PosCTL", sample_id), "PosCTL", "NegCTL")) %>%
#'     group_by(file) %>%
#'     mutate(control_group = sub(control_type, "", sample_id)) %>%
#'     mutate(control_group = sub("^_", "", control_group)) %>%
#'     mutate(control_group = sub("_$", "", control_group)) %>%
#'     mutate(sample_id = paste(control_group, control_type, sep = "_"))
#'   
#'   
#'   cell_df = bind_rows(
#'     cell_df.by_type$cohort,
#'     cell_df.by_type$control
#'   )
#' }

# library(magrittr)
# library(ggplot2)
# library(tidyverse)
# 
# .load_csv = function(file)    tryCatch({
#       dt = readr::read_csv(file)
#           }, error = function(e){
#       message("Error loading file ", file, ": ", e$message)
#       return(NULL)
#     })
# 
# .load_csv_list = function(files){
#   all_dt_l = list()
#   for(name in names(files)){
#     message("Loading data for ", name)
#     file = files[[name]]
#     stopifnot(file.exists(file))
#     dt = .load_csv(file)
#     all_dt_l[[name]] = dt
#   }
#   all_dt_l
# }
# 
# get_original_cell_data_dir = function(){
#   win_dir = r"(C:\Users\boydj\OneDrive - UVM Larner College of Medicine\Lee, Kyra C's files - VolaricDataAndScriptsForJoe\)"
#   lin_dir = "/gpfs1/home/j/r/jrboyd/VolaricDataAndScriptsForJoe/"
#   data_dir = win_dir
#   if(!dir.exists(data_dir)){
#     data_dir = lin_dir
#   }
#   stopifnot(dir.exists(data_dir))
#   data_dir
# }
# 
# load_phenocycler_summary_files = function(){
#   data_dir = get_original_cell_data_dir()
#   if(!dir.exists(data_dir)){
#     data_dir = lin_dir
#   }
#   stopifnot(dir.exists(data_dir))
#   res_files = list.files(data_dir, pattern = "Summary.+csv", recursive = TRUE,  full.names = TRUE)
#   res_files = as.list(res_files)
#   rootname = function(x){
#     x = sapply(strsplit(unlist(x), split = "/"), function(x){x[length(x)-1]})
#     x
#   }
#   names(res_files) = rootname(res_files)
#   
#   
#   # to_load = c("EBV_Pos", "EBV_PosCompanion", "EBVNeg")
#   to_load = names(res_files)
#   stopifnot(all(to_load %in% names(res_files)))
#   
#   name = "EBV_Pos"
#   all_dt_l = .load_csv_list(res_files[to_load])
#   dt = dplyr::bind_rows(all_dt_l, .id = "source") %>%
#     dplyr::mutate(Sample = sub("\\..+", "", `Image Tag`)) %>%
#     dplyr::mutate(Sample = sub("_Scan.+", "", Sample)) %>%
#     dplyr::mutate(Sample = gsub("-", "", Sample))
#   dt
# }
# 
# load_phenocycler_summary_files()
# 
# load_rnascope_summary_files = function(){
#   data_dir = get_original_cell_data_dir()
#   res_files = list.files(data_dir, pattern = "RNA.+csv", recursive = TRUE, full.names = TRUE)
#   #pivot wider results
#   res_files = res_files[!grepl("Wide", res_files)]
#   #per cell results
#   res_files = res_files[!grepl("Object", res_files)]
#   #just select single
#   res_files = res_files[grepl("RNAScopeIF_Coexpression_2026-02-11", res_files) | grepl("RNAScope_Coexpression_2026-01-12.csv", res_files)]
#   file_rename = c(
#     "RNAScope_Coexpression_2026-01-12.csv" = "4plex",
#     "RNAScopeIF_Coexpression_2026-02-11.csv" = "3plex+IF"
#   )
#   names(res_files) = file_rename[basename(res_files)]
#   res_files = as.list(res_files)
#   
#   
#   # to_load = c("EBV_Pos", "EBV_PosCompanion", "EBVNeg")
#   to_load = names(res_files)
#   stopifnot(all(to_load %in% names(res_files)))
#   
#   name = "EBV_Pos"
#   
#   all_dt_l = .load_csv_list(res_files[to_load])
#   dt = dplyr::bind_rows(all_dt_l, .id = "source") 
#   dt
# }
# 
# load_rnascope_summary_files()
# 
# 
# # all(is.na(rscope_dt$Sample) | is.na(rscope_dt$SampleNumber))
# # all(!is.na(rscope_dt$Sample) | !is.na(rscope_dt$SampleNumber))
# # rscope_dt[is.na(Sample), Sample := SampleNumber]
# # table(rscope_dt$Sample)
# # by_source = split(rscope_dt$Sample, rscope_dt$source)
# # seqsetvis::ssvFeatureUpset(by_source)
# # 
# # tmp = split(dt, dt$source)
# # tmp$`4PlexRNAScopeCellPelletSingleExpression_2026-01-12.csv`
# # tmp$`RNAScope_Coexpression_2026-01-12.csv`
# # tmp$`RNAScopeIF_Coexpression_2026-02-11.csv`
# 
# load_meta_data = function(){
#   data_dir = get_original_cell_data_dir()
#   meta_file = list.files(data_dir, pattern = "EBVStatus.xlsx", recursive = TRUE, full.names = TRUE)
#   meta_df = openxlsx::read.xlsx(meta_file)
#   meta_df$SampleRaw = meta_df$Sample
#   
#   meta_df = meta_df %>% mutate(Sample = sub("\\..+", "", SampleRaw))
#   meta_df = meta_df %>% mutate(Sample = gsub("-", "_", Sample))
#   meta_df = meta_df %>% mutate(Sample = sub("CTEBV1", "CTEBV_1", Sample))
#   meta_df = meta_df %>% mutate(Sample = sub("CTEBV2", "CTEBV_2", Sample))
#   meta_df = meta_df %>% mutate(Sample = sub("CTEBV6", "CTEBV_6", Sample))
#   meta_df = meta_df %>% mutate(Sample = sub("CTEBV9", "CTEBV_9", Sample))
#   
#   meta_df$sample = meta_df$Sample
#   
#   meta_df$SampleStripped = gsub("_", "", meta_df$Sample)
#   meta_df
# }
# 
# load_meta_data()
# 
# harmonize_phenocycler_summary_files = function(){
#   pcycler_dt = load_phenocycler_summary_files()  
#   meta_df = load_meta_data()
#   s2s = meta_df$sample
#   names(s2s) = meta_df$SampleStripped
#   
#   pcycler_dt$SampleID = s2s[pcycler_dt$Sample]
#   
#   anno_df = meta_df %>% select(SampleID = sample, EBV_Status)
#   
#   pcycler_dt = merge(pcycler_dt, anno_df, all.x = TRUE)
#   pcycler_dt = pcycler_dt %>% mutate(EBV_Status = ifelse(is.na(EBV_Status), "need info", EBV_Status))
#   pcycler_dt[]
# }
# 
# harmonize_phenocycler_summary_files()
# 
# 
# harmonize_rnascope_summary_files = function(){
#   rscope_dt = load_rnascope_summary_files()
#   
#   meta_df = load_meta_data()
#   s2s = meta_df$sample
#   names(s2s) = meta_df$SampleStripped
#   
#   rscope_dt$SampleID = s2s[rscope_dt$SampleNumber]
#   
#   
#   anno_df = meta_df %>% select(SampleID = sample, EBV_Status)
#   
#   rscope_dt = merge(rscope_dt, anno_df, all.x = TRUE)
#   rscope_dt = rscope_dt %>% mutate(EBV_Status = ifelse(is.na(EBV_Status), "need info", EBV_Status))
#   
#   rscope_dt[]
# }
# 
# harmonize_rnascope_summary_files()
# 
# # meta_full_df = openxlsx::read.xlsx("eber_status.xlsx")
# # meta_img_df = openxlsx::read.xlsx("P1_phenocycler_and_rnascope/VolaricDataAndScriptsForJoe/EBVStatus.xlsx")
# # head(meta_full_df)
# # head(meta_img_df)
# 
# # meta_img_df = meta_img_df %>% mutate(sample = sub("\\..+", "", Sample), img_status = EBV_Status, .keep = "none")
# # colnames(meta_full_df) = c("sample", "full_status")
# 
# .get_status_file = function(){
#   def_file = "/gpfs1/pi/avolaric/files_jrboyd/EBVhelpR/inst/extdata/eber_status.xlsx"
#   if(file.exists(def_file)){
#     return(def_file)
#   }
#   system.file("extdata", "eber_status.xlsx", package = "EBVhelpR")
# }
# 
# load_full_meta = function(){
#   #Updated WGS clinical annotation from Ashley
#   meta_df = openxlsx::read.xlsx(.get_status_file())
#   # meta_df = openxlsx::read.xlsx("P1_phenocycler_and_rnascope/VolaricDataAndScriptsForJoe/EBVStatus.xlsx")
#   colnames(meta_df) = c("sample", "EBER_status")
#   head(meta_df)
#   meta_df$sample = gsub("-", "_", meta_df$sample)
#   meta_df$sample %>% cbind
#   
#   meta_df$sample = gsub("[()]", "", meta_df$sample)
#   
#   #cleanup complicated sample names
#   prev_sample = grepl("previously", meta_df$sample)
#   prev_ids = strsplit(meta_df$sample[prev_sample], " previously ")[[1]]
#   prev_status = meta_df$EBER_status[prev_sample]
#   
#   
#   meta_df = rbind(
#     meta_df[!prev_sample,],
#     data.frame(sample = prev_ids, EBER_status = prev_status)
#   )
# 
#   split_ids = grepl("/", meta_df$sample)
#   split_df = meta_df[split_ids,]
#   left_df = split_df
#   right_df = split_df
#   sp = strsplit(split_df$sample, "/")
#   left_df$sample = sapply(sp, function(x)x[1])
#   right_df$sample = sapply(sp, function(x)x[2])
#   
#   meta_df = rbind(
#     meta_df[!split_ids,],
#     left_df,
#     right_df
#   )
#   
#   meta_df
# }
# 
# meta_img_df = load_meta_data()
# meta_full_df = load_full_meta()
# 
# merge(meta_img_df, meta_full_df, by = "sample", all = TRUE)
