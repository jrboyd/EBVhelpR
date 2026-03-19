.group_image_files <- function(files) {
  if(grepl("PhenocyclerReAnalysis_102025", files)){
    group = "Phenocycler"
  }else{
    group = paste0("RNAScope_", ifelse(grepl("RNAScopeIF", files), "3plex+IF", "4plex"))  
  }
  message("group is ", group)
  group
  
}


#' Title
#'
#' @returns
#' @export
#'
#' @examples
write_all_package_data = function(){
  data_dir = get_image_data_dir()
  
  rscop_object_files = dir(data_dir, pattern = "ObjectData_Clean", full.names = TRUE)
  rscop_object_files = rscop_object_files[!grepl("Sara", rscop_object_files)]
  pheno_object_files = dir(file.path(data_dir, "PhenocyclerReAnalysis_102025"), pattern = "Total_Object_Results.+csv$", recursive = TRUE, full.names = TRUE)
  
  for(file in rscop_object_files){
    write_package_data_for_file(file)
  }
  
  for(file in pheno_object_files){
    write_package_data_for_file(file)  
  }  
  
}



#' Write per-sample RNAscope cell object data
#'
#' Reads a combined RNAscope cell object CSV file and splits it into one CSV
#' per image sample, writing them to the package data directory alongside a
#' metadata summary file. Run once per input file; exits early if the metadata
#' output already exists.
#'
#' @param file Path to a combined RNAscope ObjectData CSV file.
#'
#' @return Invisibly returns the path to the metadata file written.
#' @export
write_package_data_for_file <- function(file) {
  f_group <- .group_image_files(file)
  out_dir <- file.path(.get_pkg_data_dir(), f_group)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  out_meta_file <- file.path(out_dir, paste0("metadata_", basename(file)))
  if (file.exists(out_meta_file)) {
    message("metadata file already exists, skipping: ", out_meta_file)
    return(invisible(out_meta_file))
  }
  
  message("reading input file ", file)
  obj_dat <- read.csv(file)
  
  single_tons <- vapply(obj_dat, function(x) length(unique(x)) == 1L, logical(1))
  obj_singletons <- unique(obj_dat[, single_tons, drop = FALSE])
  obj_dat <- obj_dat[, !single_tons, drop = FALSE]
  
  possible_image_loc_vars = c(
    "ImageLocation",
    "Image.Location"
  )
  loc_var_found = FALSE
  for(loc_var in possible_image_loc_vars){
    if(loc_var %in% colnames(obj_dat)){
      obj_dat$ImageLocation = obj_dat[[loc_var]]
      loc_var_found = TRUE
    }
  }
  if(!loc_var_found){
    message(paste(colnames(obj_dat), collapse = "\n"))
    stop("Could not determine image location variable")
  }
  
  obj_dat.by_image <- split(obj_dat, obj_dat$ImageLocation)
  obj_images <- names(obj_dat.by_image)
  
  if (basename(file) == "RNAScope_CellPellet_ObjectData_Cleaned_2026-01-12.csv") {
    str_extract <- function(x, delim = "_") {
      vapply(strsplit(x, delim), function(parts) {
        paste(parts[length(parts)], parts[length(parts) - 2], sep = "_")
      }, character(1))
    }
  } else {
    str_extract <- function(x) x
  }
  
  obj_names <- gsub("\\\\", "/", obj_images)
  obj_names <- basename(obj_names)
  obj_names <- sub("\\..+", "", obj_names)
  obj_names <- str_extract(obj_names)
  
  img_meta_df <- data.frame(image_name = obj_names, image = obj_images, stringsAsFactors = FALSE)
  for (col in colnames(obj_singletons)) {
    img_meta_df[[col]] <- obj_singletons[[col]]
  }
  img_meta_df[["source_file"]] <- file
  
  names(obj_names) <- obj_images
  names(obj_dat.by_image) <- obj_names[names(obj_dat.by_image)]
  
  for (name in img_meta_df$image_name) {
    message("writing ", name)
    obj_dat.sel <- obj_dat.by_image[[name]]
    obj_dat.sel$ImageLocation <- NULL
    obj_dat.sel$image_name <- name
    write.csv(obj_dat.sel, file.path(out_dir, paste0(name, ".cell_data.csv")))
  }
  
  write.csv(img_meta_df, out_meta_file)
  invisible(out_meta_file)
}