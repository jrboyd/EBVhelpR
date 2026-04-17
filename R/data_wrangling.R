.group_image_files <- function(files) {
  if(grepl("PhenocyclerReAnalysis_102025", files)){
    group = "Phenocycler"
  }else{
    group = paste0("RNAScope_", ifelse(grepl("RNAScopeIF", files), "3plex+IF", "4plex"))
  }
  message("group is ", group)
  group

}


#' Write All Package Cell Data Files
#'
#' Discovers supported RNAscope and Phenocycler object-level CSV files and
#' writes per-sample package data outputs for each file via
#' [write_package_data_for_file()].
#'
#' @return Invisibly returns `NULL` after attempting to write all package data files.
#' @examples
#' \dontrun{
#' write_all_package_data()
#' }
#' @export
write_all_package_data = function(){
  data_dir = get_original_cell_data_dir()

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
#' @examples
#' \dontrun{
#' write_package_data_for_file("/path/to/ObjectData_Clean.csv")
#' }
#' @importFrom readr read_csv
#' @export
write_package_data_for_file <- function(file) {
  f_group <- .group_image_files(file)
  out_dir <- file.path(get_wrangled_cell_data_dir(), f_group)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  out_meta_file <- file.path(out_dir, paste0("metadata_", basename(file)))
  # if (file.exists(out_meta_file)) {
  #   message("metadata file already exists, skipping: ", out_meta_file)
  #   return(invisible(out_meta_file))
  # }

  message("reading input file ", file)
  obj_dat <- readr::read_csv(file)

  single_tons <- vapply(obj_dat, function(x) length(unique(x)) == 1L, logical(1))
  obj_singletons <- unique(obj_dat[, single_tons, drop = FALSE])
  obj_dat <- obj_dat[, !single_tons, drop = FALSE]

  possible_image_loc_vars = c(
    "ImageLocation",
    "Image.Location",
    "Image Location"
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
      out_f = file.path(out_dir, paste0(name, ".cell_data.csv"))
      if(file.exists(out_f)){
          message( name, " skipped")
          next
      }
    message("writing ", name)
    obj_dat.sel <- obj_dat.by_image[[name]]
    obj_dat.sel$ImageLocation <- NULL
    obj_dat.sel$image_name <- name
    readr::write_csv(obj_dat.sel, out_f)
  }

  write.csv(img_meta_df, out_meta_file)
  invisible(out_meta_file)
}

.get_example_images_df = function(){
    f = system.file(package = "EBVhelpR", "extdata/images_on_vacc.csv", mustWork = TRUE)
    tiff_df = read.csv(f)
    tiff_df$name = NULL
    tiff_df$tiff_location = "EXAMPLE_ONLY"
    tiff_df
}


#' Title
#'
#' @returns
#' @export
#'
#' @examples
get_tiff_file_path_df = function(){
    win_dir = "Z:/FUSION DATA/AshleyVolaric"
    win_dir2 = "C:/Users/boydj/project_data/EBV_image_files"
    lin_dir = "/netfiles/volaric_research/DLBCL_EBV_detection/image_files"
    tiff_dir = win_dir
    if(!dir.exists(tiff_dir)){
        tiff_dir = lin_dir
    }
    if(!dir.exists(tiff_dir)){
        tiff_dir = win_dir2
    }
    if(!dir.exists(tiff_dir)){
        warning("Could not locate TIFF root directory. Returning example tiffs with fake paths.")
        tiff_df = .get_example_images_df()
        tiff_df$project_name = tiff_df$assay
        tiff_df$assay = project_name_to_assay[tiff_df$project_name]
    }else{
        dir_names = intersect(.get_valid_project_names(), dir(tiff_dir))
        names(dir_names) = dir_names
        tiff_files.by_project = lapply(dir_names, function(d){
            files = dir(file.path(tiff_dir, d), recursive = TRUE, pattern = "tiff?$", full.names = TRUE)
            data.frame(tiff_file = files)
        })
        tiff_df = dplyr::bind_rows(tiff_files.by_project, .id = "assay")

        tiff_df = tiff_df %>% dplyr::mutate(name = basename(tiff_file))
        tiff_df = tiff_df %>%
            dplyr::mutate(name = sub("\\..+", "", name)) %>%
            dplyr::mutate(name = sub("^[0-9]{3}_", "", name)) %>%
            dplyr::mutate(name = sub("_ ?croo?p?ped$", "", name)) %>%
            dplyr::mutate(name = gsub("-", "_", name))

        tiff_df = tiff_df %>%
            dplyr::mutate(name = ifelse(
                grepl("D_EB.+D_EB", name),
                sub("_D", " D", name),
                name
            ))
        tiff_df = tiff_df %>%
            dplyr::mutate(name = ifelse(
                grepl("CTEBV.+CTEBV", name),
                sub("_CTEBV", " CTEBV", name),
                name
            ))

        tiff_df = tiff_df %>% dplyr::group_by(assay, tiff_file) %>% dplyr::reframe(name = strsplit(name, " ")[[1]])
        tiff_df = tiff_df %>%
            dplyr::mutate(name = ifelse(
                grepl("CTEBV[0-9]", name),
                sub("CTEBV", "CTEBV_", name),
                name
            ))
        .cp_name = function(x){
            x = sub("_2$", "", x)
            x = sub("Chiung", "Chung", x)
            cells = sapply(strsplit(x, "_"), function(x)x[length(x)])
            is_control = grepl("NegCTL", x) | grepl("Control", x)
            x[is_control]
            x[!is_control]
            ifelse(is_control, paste(cells, "NegCTL", sep = '_'), cells)

        }
        tiff_df = tiff_df %>%
            dplyr::mutate(name = ifelse(
                grepl("CellPellet", name),
                .cp_name(name),
                name
            ))

        tiff_df = tiff_df %>%
            dplyr::mutate(name = sub("2468_", "", name)) %>%
            dplyr::mutate(name = sub("_Rescanned_Cropped", "", name))

        tiff_df$sample_id = tiff_df$name
        tiff_df$name = NULL

        tiff_df$project_name = tiff_df$assay
        tiff_df$assay = project_name_to_assay[tiff_df$project_name]
    }

    tiff_df$probe_control = ""
    tiff_df = tiff_df %>% dplyr::mutate(probe_control = ifelse(grepl("[Nn]eg", sample_id), "negative_probe", probe_control))
    tiff_df = tiff_df %>% dplyr::mutate(probe_control = ifelse(grepl("[Pp]os", sample_id), "positive_probe", probe_control))

    tiff_df$sample_id <- sub("_?NegCTL", "", tiff_df$sample_id)
    tiff_df$sample_id <- sub("_?PosCTL", "", tiff_df$sample_id)
    tiff_df$sample_id = sub("^JRB_", "", tiff_df$sample_id)

    #remove problematic image
    if(any(grepl("2468_D-EB-3", tiff_df$tiff_file))){
        bad_df = tiff_df %>% subset(grepl("2468_D-EB-3", tiff_file))
        warning("removing known problematic tiff : ", paste(bad_df$tiff_file, collapse = ", "))
        tiff_df = tiff_df %>% subset(!grepl("2468_D-EB-3", tiff_file))
    }


    tiff_df
}
