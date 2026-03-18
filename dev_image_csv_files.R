data_dir = "~/VolaricDataAndScriptsForJoe/"
csv_files =dir(data_dir, pattern = "RNAScope_.+csv", full.names = TRUE)
is_obj = grepl("ObjectData", csv_files)
csv_viles.is_obj = split(csv_files, is_obj)

