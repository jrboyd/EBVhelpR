library(tidyverse)
library(RBioFormats)
library(EBVhelpR)

tiff_dirs = dir("/netfiles/volaric_research/DLBCL_EBV_detection/image_files/", full.names = TRUE)
all_tiffs = dir("/netfiles/volaric_research/DLBCL_EBV_detection/image_files/", recursive = TRUE, full.names = TRUE)

df = data.frame(tiff_dir = tiff_dirs)
df = subset(df, dir.exists(tiff_dir))
df = df %>% mutate(assay = basename(tiff_dir))



tiff_df = df %>% group_by(tiff_dir, assay) %>% reframe(file = dir(tiff_dir, full.names = TRUE))


EBVhelpR::read_tiff_meta_data(tiff_df$file[1], only_series1 = FALSE)
EBVhelpR::read_tiff_meta_data(tiff_df$file[2], only_series1 = FALSE)

tmp = RBioFormats::read.metadata(tiff_df$file[1])
tmp = RBioFormats::read.metadata(tiff_df$file[2])
length(tmp)
tmp@.Data[[1]]$coreMetadata %>% length
tmp@.Data[[1]]$coreMetadata
tmp@.Data[[1]]$globalMetadata
tmp@.Data[[1]]$seriesMetadata
tmp@.Data[[1]]
tmp$seriesMetadata

f = tiff_df$file[1]
fetch_full_tiff = function(f){
  tiff_max = EBVhelpR::read_tiff_meta_data(f) %>% subset(resolutionLevel == 1)
  EBVhelpR::fetchTiffData(f, TiffRect(1, tiff_max$sizeX, 1, tiff_max$sizeY))
}

fetch_full_tiff(tiff_df$file[1])
fetch_full_tiff(tiff_df$file[2])

tiff_files = tiff_df$file
tiff_df = tiff_df %>% mutate(sample_id = sub("_cro+p.+", "", basename(file)))
names(tiff_files) = tiff_df$sample_id

f = "/netfiles/volaric_research/DLBCL_EBV_detection/image_files//RNAScopeRound1/D-EB-33.tif"

my_load_meta_data = function(f){
  message(f)
  raw_meta = tryCatch({
    RBioFormats::read.metadata(f)
  }, error = function(e){
    NULL
  })
  table_meta = tryCatch({
    EBVhelpR::read_tiff_meta_data(f)
  }, error = function(e){
    NULL
  })
  list(raw = raw_meta, table = table_meta, 
       file = f)
}

all_meta = lapply(tiff_files, my_load_meta_data)

meta_files = sapply(all_meta, function(x){
  x$file
})

can_load_meta = sapply(all_meta, function(x){
  !is.null(x$raw)
})
file.size(meta_files[which(!can_load_meta)]) %>% format(unit = "GB")


n_channels = sapply(all_meta, function(x){
  x$table$sizeC %>% max
})

f = "/netfiles/volaric_research/DLBCL_EBV_detection/image_files//Phenocycler/CellPelletSlide_Control_Scan1_Phenocycler_GM12878.ome.tiff"
tmp = RBioFormats::read.metadata(f)
length(tmp[[1]]$coreMetadata)
# tmp[[1]]$coreMetadata$


img_data = read.image(f, resolution = 1)
dim(img_data)
plot(0:1, 0:1)
i = 7
rasterImage(img_data[,,i] / max(img_data[,,i]), 0, 0, 1, 1)


tmp[[2]]

all_meta$`201_D-EB-2`

all_meta[[20]]$coreMetadata


try_load_meta = function(f){
  tryCatch({
    EBVhelpR::read_tiff_meta_data(f, only_series1 = FALSE)  
  }, error = function(e){
    return(data.frame(series = NA, sizeX = NA, sizeY = NA, sizeC = NA, resolutionLevel = NA))
  })
}

try_load_data = function(f, xmax, ymax){
  tryCatch({
    tdat = EBVhelpR::fetchTiffData(f, TiffRect(1, xmax, 1, ymax))  
    return(TRUE)
  }, error = function(e){
    return(FALSE)
  })
}

meta_df = tiff_df %>% group_by(tiff_dir, assay, file) %>% reframe(try_load_meta(file))
meta_df = meta_df %>% group_by(tiff_dir, assay, file) %>% summarise(n_series = length(unique(series)), sizeX = max(sizeX), sizeY = max(sizeY), n_resolutions = max(resolutionLevel), sizeC = max(sizeC))
meta_df = meta_df %>% group_by(file) %>% mutate(can_load = try_load_data(file, sizeX, sizeY))
meta_df

#### load object data ####
library(data.table)
objdata_files = dir("~/VolaricDataAndScriptsForJoe/", pattern = "ObjectData", full.names = TRUE)
names(objdata_files) = basename(objdata_files)
objdata_files = objdata_files[!grepl("Sara", objdata_files)]
all_objdata = lapply(objdata_files, function(f){
  dt = tryCatch({
    fread(f)
  }, error = function(e){
    NULL
  })
  dt
})

all_objdata = all_objdata[!sapply(all_objdata, is.null)]
names(all_objdata)
all_objdata.l = lapply(all_objdata, function(x){
  split(x, x[[1]])
})

x = all_objdata.l[[1]]
all_objdata.l = lapply(all_objdata.l, function(x){
  xnames = sapply(strsplit(names(x) , "\\\\"), function(y)y[length(y)])
  names(x) = sub("\\..+", "", xnames)
  x
})

sapply(all_objdata.l, names)

# focus on scopeIF
table(tiff_df$assay)
sel_tiff_df = subset(tiff_df, assay == "RNAScopeRound2")
sel_tiff_df$sample_id = sel_tiff_df$sample_id %>% sub("\\..+", "", .)
sel_tiff_df$sample_id = sel_tiff_df$sample_id %>% sub("_cro.+", "", .)
intersect(sel_tiff_df$sample_id, names(all_objdata.l$RNAScopeIF_ObjectData.csv))
names(all_objdata.l$RNAScopeIF_ObjectData.csv)

setdiff(sel_tiff_df$sample_id, names(all_objdata.l$RNAScopeIF_ObjectData.csv))

valid_ids = intersect(sel_tiff_df$sample_id, names(all_objdata.l$RNAScopeIF_ObjectData.csv))
cell_data = all_objdata.l$RNAScopeIF_ObjectData.csv[valid_ids]  
sel_tiff_df = subset(sel_tiff_df, sample_id %in% valid_ids)

f = sel_tiff_df$file[1]
sel_id = subset(sel_tiff_df, file == f)$sample_id
sel_data = cell_data[[sel_id]]
dim(sel_data)

full_image = fetch_full_tiff(f)
full_rect = full_image@rect
sel_rect = full_rect
select_cells = function(sel_rect){
  cell_data = subset(sel_data, XMax > sel_rect@xmin & XMin < sel_rect@xmax & YMax > sel_rect@ymin & YMin < sel_rect@ymax)
  cell_data %>% select(`Object Id`, XMin, XMax, YMin, YMax)
}
plot_cells = select_cells(full_rect)
sel_data
plot(full_image) + annotate("point", x = (plot_cells$XMin + plot_cells$XMax)/2, y = (plot_cells$YMin + plot_cells$YMax) / 2, color = "green")

my_zoom = function(tiff_dat, zoom_rect){
  p = tiff_dat@plots$channels
  p_anno = p %>% rect_annotate(zoom_rect)
  zoom_dat = fetchTiffData(tiff_dat$tiff_path, rect = zoom_rect, precalc_max = tiff_dat@precalc_max)
  zoom_dat@plots[["original"]] = p_anno
  zoom_dat@plots[["assembly"]] = cowplot::plot_grid(p_anno, zoom_dat@plots$channels)
  zoom_dat@activePlot = "assembly"
  zoom_dat
}

#### pheno tiffs ####
pc_df = subset(tiff_df, assay == "Phenocycler")
pc_df.meta = lapply(pc_df$file, function(f){
  df = try_load_meta(f)
  if(!is.null(df)){
    df = subset(df, resolutionLevel == 1)
  }
  df
})
pc_df = cbind(pc_df, do.call(rbind, pc_df.meta))
res = system(paste0("md5sum '", f, "'"), intern = TRUE)
mres = strsplit(res, " .+")[[1]][1]
md5sums = sapply(pc_df$file, function(f){
  res = system(paste("md5sum", f), intern = TRUE)
  mres = strsplit(res, " .+")[[1]][1]
  mres
})
tmp = read.image(f, resolution = 4, filter.metadata = FALSE)
tmp@metadata$coreMetadata


