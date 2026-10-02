test_that("the registry is internally consistent", {
  locs <- ebv_data_locations()
  src  <- ebv_data_sources()
  expect_false(anyNA(locs$priority))
  # one priority, role and mirror_subdir per location
  per <- unique(locs[, c("location_id", "priority", "role", "mirror_subdir")])
  expect_equal(anyDuplicated(per$location_id), 0)
  expect_equal(anyDuplicated(src$dataset_id), 0)
  expect_true(all(src$kind %in% c("file", "dir")))
  expect_true(all(src$tier %in% c("mirror", "split", "reference")))
  ids <- unique(locs$location_id)
  has_copy <- rowSums(!is.na(as.matrix(src[, ids]))) > 0
  expect_true(all(has_copy), info = paste(src$dataset_id[!has_copy], collapse = ", "))
  # relative paths only, forward slashes
  rels <- unlist(src[, ids]); rels <- rels[!is.na(rels)]
  expect_false(any(grepl("\\\\", rels)))
  expect_false(any(grepl("^(/|[A-Za-z]:)", rels)))
})

# A fake world: one source location under a temp dir, every other location
# pointed at a path that does not exist, and a temp mirror.
local_fake_world <- function(env = parent.frame()) {
  root   <- withr::local_tempdir(.local_envir = env)
  mirror <- withr::local_tempdir(.local_envir = env)
  nowhere <- file.path(root, "does-not-exist")
  ids <- unique(ebv_data_locations()$location_id)
  vars <- stats::setNames(rep(nowhere, length(ids)), paste0("EBVHELPER_ROOT_", toupper(ids)))
  vars[["EBVHELPER_ROOT_LABSHARE"]] <- root
  vars[["EBVHELPER_MIRROR_DIR"]] <- mirror
  vars[["EBVHELPER_CLEANUP_LOG"]] <- ""
  withr::local_envvar(vars, .local_envir = env)
  list(root = root, mirror = mirror)
}

test_that("ebv_location_root honours the environment override", {
  w <- local_fake_world()
  expect_equal(ebv_location_root("labshare"), w$root)
})

test_that("ebv_path resolves a source, then the mirror after fetching, and fails loudly", {
  w <- local_fake_world()
  rel <- ebv_data_sources()
  rel <- rel$labshare[rel$dataset_id == "ihc_clean_2025_09_03"]
  f <- file.path(w$root, rel)
  dir.create(dirname(f), recursive = TRUE)
  writeLines(c("Sample,EBNA2", "DEB1,0"), f)

  # the override is checked before the real candidates, so this is the fake root
  p <- ebv_path("ihc_clean_2025_09_03", use_mirror = FALSE)
  expect_equal(normalizePath(p), normalizePath(f))

  res <- ebv_mirror_fetch("ihc_clean_2025_09_03", hash = TRUE)
  expect_equal(res$result, "copied")
  mp <- ebv_path("ihc_clean_2025_09_03")
  expect_true(startsWith(normalizePath(mp), normalizePath(w$mirror)))
  expect_equal(readLines(mp), readLines(f))
  man <- ebv_mirror_manifest()
  expect_equal(nrow(man), 1)
  expect_equal(man$md5, unname(tools::md5sum(mp)))

  # unchanged on a second fetch; a changed source is NOT overwritten, only logged
  expect_equal(ebv_mirror_fetch("ihc_clean_2025_09_03")$result, "unchanged")
  writeLines(c("Sample,EBNA2", "DEB1,0", "DEB2,90"), f)
  Sys.setFileTime(f, Sys.time() + 60)
  expect_equal(ebv_mirror_fetch("ihc_clean_2025_09_03")$result, "differs")
  expect_equal(length(readLines(mp)), 2)
  log <- utils::read.delim(ebv_cleanup_log(), stringsAsFactors = FALSE)
  expect_equal(nrow(log), 1)
  expect_match(log$reason, "not overwritten")

  # gone everywhere -> informative error
  unlink(f)
  unlink(mp)
  expect_error(ebv_path("ihc_clean_2025_09_03"), "not available on this machine")
  expect_true(is.na(ebv_path("ihc_clean_2025_09_03", must_work = FALSE)))
})

test_that("reference and never-mirror datasets are refused", {
  w <- local_fake_world()
  expect_error(ebv_mirror_fetch("img_tree"), "reference")
  expect_error(ebv_mirror_fetch("ihc_master_case_list", allow_reference = TRUE), "Never mirror")
})

test_that("ebv_log_cleanup appends once per path and reason, and touches nothing", {
  w <- local_fake_world()
  target <- file.path(w$root, "dup.csv"); writeLines("x", target)
  expect_true(ebv_log_cleanup(target, "duplicate", keep = "elsewhere", location = "labshare"))
  expect_false(ebv_log_cleanup(target, "duplicate"))
  expect_true(ebv_log_cleanup(target, "truncated"))
  log <- utils::read.delim(ebv_cleanup_log(), stringsAsFactors = FALSE)
  expect_equal(nrow(log), 2)
  expect_equal(log$bytes[1], file.size(target))
  expect_true(file.exists(target))
})

test_that("a NUL-filled mirror copy is not trusted", {
  w <- local_fake_world()
  rel <- ebv_data_sources(); rel <- rel$labshare[rel$dataset_id == "ihc_long_2025_09_04"]
  mp <- file.path(w$mirror, "cache", "files.med.uvm.edu/shared/Centers/CBSR/PI/Volaric", rel)
  dir.create(dirname(mp), recursive = TRUE)
  writeBin(raw(1024), mp)
  expect_true(is.na(ebv_path("ihc_long_2025_09_04", must_work = FALSE)) ||
              !startsWith(normalizePath(ebv_path("ihc_long_2025_09_04")), normalizePath(w$mirror)))
})
