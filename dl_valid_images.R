img_file_df = read.csv("inst/extdata/images_on_vacc.csv")

valid_cell_pngs = dir("validated_cells/", recursive = TRUE, pattern = ".png$", full.names = TRUE)
valid_df = data.frame(valid_file = valid_cell_pngs)
valid_df = valid_df %>% mutate(assay = basename(dirname(valid_file)) %>% sub("_JRB", "", .), unique_id = basename(valid_file) %>% sub("00_validate_", "", .) %>% sub(".png", "", .))
valid_ids.by_assay = split(valid_df$unique_id, valid_df$assay)


sel_df = valid_df %>% select(sample_id = unique_id, assay)
sel_df$assay = EBVhelpR:::assay_to_project_name[sel_df$assay]

sel_df = merge(img_file_df, sel_df)
sel_df = sel_df %>% mutate(fetch_name = file.path(basename(dirname(tiff_file)), basename(tiff_file)))
cbind(sel_df$fetch_name)
message(paste(sel_df$fetch_name, collapse = "\n"))
#append to /cygdrive/c/Users/boydj/project_data/EBV_image_files/images.txt
