# EBVhelpR — known issues

Compiled 2026-09-09, updated 2026-09-15, from a QC pass over the EBV/DLBCL
Phenocycler, RNAscope and WGS data (`paper1_enhanced_EBV_detection`). Every item
was verified against the source at the line given. Ordered by what it costs.

---

## Correctness — these have changed published numbers

### 1. `set_selected_unique_ids()` does not update `selected_sample_ids`

`R/query_class.R:333` sets only `selected_unique_ids`, but
`get_query_summary_df()` (`R/query_class.R:259`) filters on **both**
`selected_sample_ids` and `selected_unique_ids`. A sample present in
`summary_df` but absent from `all_cell_files_df` therefore passes the
overlap/upset selection and then silently disappears from every downstream
figure.

**Observed:** `D_EB_33` is in the 4-way assay overlap set but every Phenocycler
panel is n = 29 rather than 30. The underlying cause is legitimate — the warning
`removing D_EB_12_D_EB_33_Scan5, multiple samples in same image` — but the
overlap logic never learns about it.

**Fix:** either have the setter derive `selected_sample_ids` from the supplied
unique ids, or reduce `get_query_summary_df()` to filtering on one key. Emit a
warning when a requested id is dropped.

### 2. `selected_unique_ids` is initialised from `sample_id`, not `unique_id`

`R/query_class.R:108`:

```r
selected_sample_ids  = all_cell_files_df$sample_id %>% unique   # :102, correct
selected_unique_ids  = all_cell_files_df$sample_id %>% unique   # :108, should be unique_id
```

Line 108 runs *after* `.df_prep()` has created `unique_id`, so the column is
available. As written, any probe-control row (`unique_id` = `"X negative_probe"`)
can never be selected, and it compounds issue 1.

### 3. `plot_wgs_pileup_heatmap()` silently collapses out-of-level samples

In `R/wgs.R`, `factor(plot_df$sample, levels = lvl)` maps any sample not present
in `wgs_count_summary` to `NA`, and `geom_tile` then draws all of them stacked at
one y position.

**Observed:** an "NA" row at the top of
`run_control/02_initial_wgs_profiles_heatmap.png` containing 37 overplotted
patient samples.

**Fix:** drop non-matching samples with a warning, or stop. Silently plotting
them as a single row is the worst option.

---

## Silent failures

### 4. `separate()` called without the `tidyr::` prefix

`R/rnascope.R:232`, inside `load_cell_source_files()`. Errors with
`could not find function "separate"` unless the caller has attached tidyverse.
It only works today because the figure scripts call `library(tidyverse)`.
Hit twice during this QC pass.

### 5. `cell_counts.txt` parse keeps the trailing `total` line

`R/rnascope.R:229-232`. The `wc -l`-style summary line becomes a data row with
`file == "total"`.

### 6. `load_meta_data()` validates two columns then assigns three names

`R/metadata.R:51` checks `ncol(meta_df) < 2`; `R/metadata.R:55` assigns
`colnames(meta_df)[1:3]`. A two-column sheet passes the check and fails
confusingly one line later.

### 7. Status-bar fills go `NA` for unrecognised statuses

`status_colors[anno_dt$EBER_status]` in the heatmap annotation draws transparent
rectangles rather than raising an error when a status has no colour.

---

## Portability

### 8. Hardcoded absolute paths with the username baked in

- `R/wgs.R:60-61` — WGS root
- `R/io_utils.R:39-40` — original cell data root
- `R/data_wrangling.R:156` — image root

Only `get_original_cell_data_dir()` offers an environment override
(`EBVHELPER_DATA_DIR`). The WGS and image roots should get the same treatment
(`EBVHELPER_WGS_DIR`, `EBVHELPER_IMAGE_DIR`). Images have already been moved to
`G:/project_data/EBV_image_files` once for disk space, which the hardcoded path
does not know about.

### 9. `.get_status_file()` prefers a cluster path over the packaged copy

`R/io_utils.R:57` checks
`/gpfs1/pi/avolaric/files_jrboyd/EBVhelpR/inst/extdata/eber_status.xlsx` first.
On the cluster the metadata can therefore differ from the installed package's
copy with no indication, so the same code gives different EBER calls in two
places.

---

## Packaging

### 10. `License: What license is it under?`

`DESCRIPTION:13`, still the template placeholder.

### 11. `TiffPlotR (>= 0.1.3)` in Imports with no `Remotes:` field

A fresh install cannot resolve it.

### 12. Uncommitted source that downstream analysis now depends on

`R/wgs.R` and `R/cell_query_workflow.R` are modified, three R files deleted, most
of `man/` modified. The installed 0.1.5 build (2026-05-12) contains these edits;
the last commit does not.

`paper1_enhanced_EBV_detection/despike_wgs.R` requires the `...` passthrough and
the `if (smooth_n > 1)` guard in `load_wgs_bigwig_pileup()`. **A reinstall from
HEAD breaks the de-spiking and everything built on it.** Worth committing first.

---

## Semantics worth documenting rather than changing

### 13. `load_wgs_bigwig_pileup()` shifts x to start at zero

`plot_dt$x <- plot_dt$x - min(plot_dt$x)` makes `x` a window-centre offset rather
than a reference coordinate — roughly half a window off, and dependent on the
fetched sample set. It matters because the GFF gene track is overlaid in true
coordinates. The offset is ~25 bp against 171 kb, so harmless in practice but
undocumented.

### 14. Technical replicates are pooled silently

`.harmonize_bigwig_sample_ids()` strips `_1`/`_2` suffixes and
`load_wgs_count_summary()` then sums reads across the collapsed id. Load-bearing
for `D_EB_65/66/67`, currently producing no collisions (51 files → 51 ids), but
nothing says so.

### 15. `EBV_OPAL_DECODE` declares Opal 690 for both RNAscope panels

Opal 690 is **not** in either per-cell export, so only three of four probes reach
the user per cell. This is deliberate — the channel was judged unusable — and the
images confirm it, with two different failures:

- **3-plex+IF (EBNA1 antibody):** Opal 690 is the channel most correlated with
  the dedicated autofluorescence channel in all three images tested (Spearman
  rho 0.33–0.52, partial on DAPI; two to ten times any other Opal), and it has no
  zero-pixel population, sitting on a floor of 2.5–2.9 where Opal 520/570 sit at
  exactly 0.00. It is measuring autofluorescence.
- **4-plex (EBNA3):** simply empty — dynamic range 1.2× above background and
  rho −0.002 against AF, while Opal 520 (EBER1) in the same image gives 42×.

**Fix:** annotate the decode table so downstream code cannot assume EBNA3 or the
EBNA1 antibody exist, rather than removing them (the channels are in the images).

### 16. Composite column names lose `+` and `−` — a correctness hazard

Upgraded from a documentation note after verification. The Phenocycler export
carries `CD4+TCells` and `CD4-Tcells`; `make.names()` maps **both `+` and `−` to
`.`**, so they collide and survive only as `CD4.TCells` and `CD4.Tcells`,
distinguished by letter case alone. Likewise `CD8+Tcells` / `CD8-Tcells` become
`CD8.Tcells` and `CD8.Tcells.1`.

Verified by assertion across all four sections tested: each pair sums exactly to
the CD3⁺ count, i.e. **they are a population and its complement.** One case typo
analyses the opposite cells. Also confirmed: `PDL1.Bcells` is PDL1 ∩ **Pax5**,
not PDL1 ∩ CD20.

**Fix:** rename on import in `load_phenocycler_summary_files()` to
`CD4_pos_Tcells` / `CD4_neg_Tcells` etc., and assert the complement relationship.

---

## Metadata

### 17. `D-EB-11` reclassified to Positive — installed copy is stale

`inst/extdata/eber_status.xlsx` was edited 2026-09-15: `D-EB-11` changed from
`Not performed` to `Positive`, with a note recorded in the Notes column. Backup
at `inst/extdata/eber_status.backup_20260915.xlsx`. The sheet now holds
61 Negative and 34 Positive, and **no `Not performed` rows at all**.

Two consequences:

- `load_meta_data()` resolves through `system.file()`, i.e. the **installed**
  package, so the change is not visible to any analysis until EBVhelpR is
  reinstalled.
- `D-EB-11` has WGS *and* all three assays (Phenocycler, RNAScope_4plex,
  RNAScope_3plex+IF). It was previously excluded everywhere by the
  `EBER_status %in% c("Positive", "Negative")` filter. On reinstall it enters the
  WGS Pos/Neg set (50 → 51) and the assay-overlap cohort (30 → 31), and it is a
  Positive case. Cached derived data — including
  `presentation_cache/despiked_wgs_coverage.Rds` — was built without it and needs
  refreshing.

Code paths handling `"Not performed"` are now unexercised but should be kept;
the value may return with new cases.

---

## Suggested order

1. Items **12** (commit the source), **1** and **2** — 12 protects the work built
   on the current build; 1 and 2 are the only items that have altered a published
   number.
2. Item **17**: reinstall, then re-run anything cached.
3. Item **16**, then **4**–**7** — cheap, and each removes a silent failure.
4. Items **8**–**11**, **13**–**15** as convenient.
