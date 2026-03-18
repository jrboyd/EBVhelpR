library(magrittr)

EBVhelpR::load_phenocycler_summary_files()
EBVhelpR::harmonize_phenocycler_summary_files()

rscop_summary = EBVhelpR::load_rnascope_summary_files()
rscop_summary$source %>% table
EBVhelpR::harmonize_rnascope_summary_files()

data_dir = EBVhelpR:::.get_data_dir()

.get_data_dir = EBVhelpR:::.get_data_dir

.get_pkg_data_dir = function(){
  main_dir = .get_data_dir()
  file.path(main_dir, "../EBVhelpR_data")
}
pkg_data_dir = .get_pkg_data_dir()

rscop_object_files = dir(data_dir, pattern = "ObjectData_Clean", full.names = TRUE)
rscop_object_files = rscop_object_files[!grepl("Sara", rscop_object_files)]
pheno_object_files = dir(file.path(data_dir, "PhenocyclerReAnalysis_102025"), pattern = "Total_Object_Results.+csv$", recursive = TRUE)

group_rscope_files = function(files){
  paste0("RNAScope_", ifelse(grepl("RNAScopeIF", files), "3plex+IF", "4plex"))
}

file = rscop_object_files[1]


write_package_data_rnascope_data = function(file){
  f_group = group_rscope_files(file)
  
  out_dir = file.path(.get_pkg_data_dir(), f_group)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  meta_f_name = paste0("metadata_", basename(file))
  out_meta_file = file.path(out_dir, meta_f_name)
  if(file.exists(out_meta_file)){
    message("prematurely exitting, meta file exists: ", out_meta_file)
    return(out_meta_file)
  }
  message("reading input file ", file)
  obj_dat= read.csv(file)
  single_tons = sapply(obj_dat, function(x){
    length(unique(x)) == 1
  })
  obj_singletons = unique(obj_dat[, single_tons])
  obj_dat = obj_dat[, !single_tons]
  obj_dat.by_image = split(obj_dat, obj_dat$ImageLocation)
  obj_images = names(obj_dat.by_image)
  if(basename(file) == "RNAScope_CellPellet_ObjectData_Cleaned_2026-01-12.csv"){
    str_extract = function(x, delim = "_"){
      sapply(strsplit(x, delim), function(x){
        paste(x[length(x)], x[length(x)-2], sep = '_')
      })
    }  
  }else if(basename(file) == "RNAScope_ObjectData_Cleaned_01122026.csv"){
    # browser()
    str_extract = function(x){x}
    # str_extract = function(x, delim = "_"){
    #   sapply(strsplit(x, delim), function(x){
    #     paste(x[length(x)], x[length(x)-2], sep = '_')
    #   })
    # }
  }else{
    # browser()
    str_extract = function(x){x}
  }
  
  obj_names = obj_images %>% gsub("\\\\", "/", .) %>% basename %>% sub("\\..+", "", .) %>% str_extract
  # normalizePath(obj_images, winslash = "\\")
  obj_names
  img_meta_df = data.frame(image_name = obj_names, image = obj_images)
  for(name in colnames(obj_singletons)){
    img_meta_df[[name]] = obj_singletons[[name]]
  }
  img_meta_df[["source_file"]] = file
  
  names(obj_names) = obj_images
  names(obj_dat.by_image) = obj_names[names(obj_dat.by_image)]
  
  name = img_meta_df$image_name[1]
  for(name in img_meta_df$image_name){
    message("writing ", name)
    obj_dat.sel = obj_dat.by_image[[name]]
    obj_dat.sel$ImageLocation = NULL
    obj_dat.sel$image_name = name

    f_name = paste0(name, ".cell_data.csv")
    out_f = file.path(out_dir, f_name)

    write.csv(obj_dat.sel, out_f)
  }
  
  write.csv(img_meta_df, out_meta_file)
}

for(file in rscop_object_files){
  write_package_data_rnascope_data(file)  
}

pkg_data_dir = .get_pkg_data_dir()

meta_df = EBVhelpR::load_meta_data()

all_cell_data_files = list.files(pkg_data_dir, pattern = ".cell_data.csv$", recursive = TRUE)

cell_meta_files = list.files(pkg_data_dir, pattern = "meta.+csv$", recursive = TRUE, full.names = TRUE)
lapply(cell_meta_files, read.csv)

cell_df = data.frame(file = all_cell_data_files)
library(tidyverse)
cell_df = cell_df %>% separate(file, sep = "/", into = c("assay", "image_name"), remove = FALSE) %>% mutate(image_name = sub("\\..+", "", image_name))
cell_df$name = cell_df$image_name
cell_df = cell_df %>% 
  mutate(name = sub("_ ?[cC]ro.+", "", name)) %>% 
  mutate(name = sub("^[0-9]+_", "", name)) %>% 
  mutate(name = sub("EBER-LMP1-EBNA1", "PosCTL", name)) %>% 
  mutate(name = sub("_Rescanned", "", name)) %>% 
  mutate(name = gsub("_Rescanned", "", name))
cell_df = cell_df %>% mutate(sample_type = "cohort") %>%
  mutate(sample_type = ifelse(grepl("CellPellet", name), "pellet", sample_type)) %>%
  mutate(sample_type = ifelse(grepl("CTL", name), "control", sample_type)) 


intersect(meta_df$sample, cell_df$name)
cell_df

dir(pkg_data_dir, full.names = TRUE) %>% dir

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
# .get_data_dir = function(){
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
#   data_dir = .get_data_dir()
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
#   data_dir = .get_data_dir()
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
#   data_dir = .get_data_dir()
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
