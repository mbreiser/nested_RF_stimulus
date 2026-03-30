# Codebase Review: nested_RF_stimulus

**Date:** 2026-02-17
**Repo:** nested_RF_stimulus (Reiser Lab, Janelia)
**Files:** ~111 MATLAB source files
**Purpose:** Two-protocol pipeline for receptive field (RF) and direction selectivity (DS) mapping in *Drosophila* T4/T5 neurons using the G4 LED arena.

---

## 1. Architecture Overview

The codebase is split into three major subsystems:

```
src/
├── stimulus_generation/    Builds patterns + position functions for G4 arena
│   ├── pattern/            Flash and bar pattern generators
│   ├── posn_fn/            Position function generators
│   └── helper/             Frame↔coordinate conversion, stimulus utilities
├── protocol_generation/    Assembles and runs Protocol 2 experiments
└── analysis/               Post-hoc data processing and visualization
    ├── protocol2/          Single-experiment and batch analysis pipelines
    │   └── pipeline/       Core parsing: parse_bar_data, parse_flash_data, etc.
    ├── plotting/           ~25 specialised plot functions
    ├── helper/             Stats, colormaps, Gaussian fitting, etc.
    ├── analyse_bar_DS/     Population-level bar DS analysis
    ├── results_analysis/   Table builders for aggregated results
    └── quality_check/      Recording quality assessment
```

**Workflow:**
Protocol 1 (pre-made) → identify peak frame → `generate_protocol2()` → run on G4 → `process_protocol2()` → figures + .mat results

---

## 2. Data Input Strategy

### 2.1 Raw Data Ingested

| Format | Source | Loaded by | Notes |
|--------|--------|-----------|-------|
| `.mat` (G4_TDMS_Log) | G4 arena recording | `load_protocol2_data.m` | Contains `Log.ADC.Volts` — row 1 = frame data, row 2 = voltage (×10 scaling) |
| `.mat` (currentExp) | Protocol generation | `process_protocol2.m`, `analyze_single_experiment.m` | Metadata: Frame, Age, Strain, Side, pattern/func order |
| `.mat` (params) | Stimulus generation | `load_protocol2_data.m` | Stimulus parameters: x, y, on_off, crop sizes |
| `.mat` (pfnparam) | Position functions | `load_protocol2_data.m` | Function parameters: frame count, duration |
| `.mat` (bar_lut) | Lookup table | `analyze_single_experiment.m` | Maps pattern IDs to directions/orientations |
| `.mat` (metadata) | Experiment metadata | `extract_metadata_from_folders.m` | Per-experiment fly info (genotype, age, sex) |
| `.mat` (bar_results_*.mat) | Previous analysis | `combine_bar_results.m` | Aggregated bar DS results |
| `.mat` (rf_results_*.mat) | Previous analysis | `make_rf_res_table.m` | Aggregated RF results |
| `.mat` (qual_res_*.mat) | Quality check | `make_qual_res_table.m` | Recording quality metrics |

### 2.2 Data Loading Patterns

- **Primary method:** `load(filename, 'varname')` — used throughout (~40+ instances)
- **Path discovery:** `dir('G4_TDMS*')`, `dir('*.mat')` — glob-based file discovery
- **Path construction:** Mix of `fullfile()` (good) and string concatenation with hardcoded separators (fragile)
- **Working directory reliance:** `load_protocol2_data.m` uses `cd()` to navigate into subfolders — side-effect-prone

---

## 3. Data Output Strategy

### 3.1 Results Saved

| Output | Format | Saved by | Contents |
|--------|--------|----------|----------|
| `bar_results_*.mat` | .mat | `process_bars_p2.m` | `bar_results` struct: DSI, FWHM, CV, vector sum, max/min voltage, PD angle |
| `rf_results_*.mat` | .mat | `process_flash_p2.m` | `rf_results` struct: Gaussian fit params (σx, σy, R²), spatial response maps, variance metrics |
| `qual_res_*.mat` | .mat | `quality_analysis.m` | Quality table: median/var/max/min voltage |
| `.pfn` binary | Binary | `save_function_G4_test.m` | G4 position function files for the arena controller |
| Google Sheets row | HTTP POST | `export_to_google_sheets.m` | Metadata via Google Form submission (Date, Time, Strain, Age, Frame, Side, Drug, Notes) |

### 3.2 Data Saving Patterns

- **Primary method:** `save(fullfile(...), 'var1', 'var2')` — standard MATLAB .mat
- **Naming convention:** `<type>_<strain>_<on_off>_<date>_<time>.mat`
- **No CSV/Excel export** for final results (all .mat); table assembly happens in `make_*_res_table.m` functions but output stays in MATLAB workspace
- **Binary I/O:** `save_function_G4_test.m` uses `fopen/fwrite/fclose` without error checking

---

## 4. Plot/Figure Output Strategy

### 4.1 Plot Types Generated

| Plot Type | Function(s) | Format |
|-----------|-------------|--------|
| Polar plots (DS tuning curves) | `plot_polar_with_arrow.m`, `plot_slow_bar_sweep_polar.m`, `plot_PD_vs_angle_polar.m`, `plot_OD_vs_angle_polar.m` | Vector PDF |
| Timeseries + embedded polar | `plot_timeseries_polar_bars.m`, `plot_timeseries_polar_bars_pharma.m` | Vector PDF |
| RF heatmaps | `plot_heatmap_flash_responses.m`, `plot_bar_flash_heatmap.m` | Vector PDF |
| Gaussian RF contours | `gaussian_RF_estimate.m` | Displayed inline |
| 1×11 bar flash profiles | `plot_bar_flash_1x11.m` | Vector PDF |
| Population polar overlays | `plot_polar_population.m`, `polar_mean_by_strain.m` | Vector PDF |
| Linear timeseries grids | `plot_timeseries_ind_cells.m`, `plot_linear_timeseries_16D.m` | PDF |
| Quality-check diagnostic plots | `plot_quality_check_*.m` (4 functions) | Displayed / optional save |
| Boxplots | `generate_boxplot_data_tbl.m` | Displayed |

### 4.2 Figure Saving Patterns

- **Primary method:** `exportgraphics(fig, path, 'ContentType', 'vector', 'BackgroundColor', 'none')` — high-quality vector PDFs
- **Secondary:** `saveas(gcf, path)` and `exportgraphics(..., 'ContentType', 'image', 'Resolution', 300)` in newer code
- **Output directories:** Created on-the-fly with `mkdir()` under `PROJECT_ROOT/figures/`
- **No PNG/SVG export** — exclusively PDF output

---

## 5. TODOs and Incomplete Work

### 5.1 Explicit TODOs

| File | Line | TODO |
|------|------|------|
| `test_protocol_length_diff_param_durations.m` | 42 | `% TODO - update to make these initial patterns from script too....` |
| `parse_flash_data.m` | 218 | `% % % % % % % TODO - - -- update this based on the 1000` (garbled, appears abandoned) |
| `plot_timeseries_polar_bars.m` | 154 | `% % % TODO - - - Add lines when the stimulus starts / stops? - 9000` |
| `plot_timeseries_polar_bars_pharma.m` | 159 | `% % % TODO - - - Add lines when the stimulus starts / stops? - 9000` (duplicate of above) |

### 5.2 Deprecated / Superseded Code

| File | Status |
|------|--------|
| `generate_protocol2_stimuli.m` | Line 7: *"This is an earlier version... recommends using GENERATE_PROTOCOL2 instead"* |
| `flash4pix_analysis.m` | Script-style file with `close all; clear;` — ~130 lines of commented-out code (lines 436–557, 621–653). Appears to be an earlier interactive analysis workflow. |
| `quality_analysis.m` | Script with `close all; clear;` — hardcoded paths. Earlier version of quality pipeline. |
| `plot_bar_responses.m` | Extensive commented-out blocks (lines 45–489); hardcoded user paths. |

### 5.3 Commented-Out Numbering Mismatch

In `generate_protocol2.m`, steps 5 and 6 are numbered twice:
- Lines 125–129: Steps 5 and 6 for bar flash generation
- Lines 131–135: Steps 5 and 6 again for protocol assembly and running

---

## 6. Potential Bugs and Issues

### 6.1 High Priority

| Issue | File:Line | Detail |
|-------|-----------|--------|
| **Hardcoded Windows path** | `process_protocol2.m:42` | `PROJECT_ROOT = "C:\matlabroot\G4_Protocols\nested_RF_protocol2"` — will fail on macOS/Linux. Same in `generate_protocol2.m:92`. |
| **Hardcoded user paths** | `plot_bar_responses.m:27–28`, `quality_analysis.m:35` | Paths contain `/Users/burnettl/...` — will fail for other users. |
| **`cd()` side effects** | `load_protocol2_data.m:33,47,57,69` | Function changes working directory 4 times via `cd()`. If an error occurs partway through, the caller's working directory is left in an inconsistent state. `analyze_single_experiment.m` mitigates this with `onCleanup`, but `process_protocol2.m` does not. |
| **File I/O without error check** | `save_function_G4_test.m:83–85` | `fopen` result not checked; if it returns -1, `fwrite` and `fclose` will error or silently fail. |
| **Silent catch blocks** | `extract_metadata_from_folders.m:158,206`, `plot_PD_vs_angle_polar.m:59`, `plot_OD_vs_angle_polar.m:56`, `combine_bar_results.m:285` | Bare `catch` without exception variable — errors are silently swallowed. Debugging will be difficult. |

### 6.2 Medium Priority

| Issue | File:Line | Detail |
|-------|-----------|--------|
| **`load_protocol2_data.m` always loads first param file** | Line 65 | `load(param_file(1).name, 'params')` — comment acknowledges this only loads 4px params. If 6px param file is listed first (OS-dependent sort), wrong params could be loaded. |
| **Hardcoded index 144 in colormap** | `inferno.m:262` | `hsv(144:end,1)=hsv(144:end,1)+1; % hardcoded` — unexplained magic number; fragile if colormap resolution changes. |
| **Path separator detection** | `load_protocol2_data.m:35–39` | Uses `contains(exp_folder, '/')` vs `'\'` to split — will fail on mixed-separator paths (e.g., `fullfile` on Windows can produce `/`). Better to use `fileparts()` or `filesep`. |
| **`idx` extraction assumes fixed structure** | `process_flash_p2.m:82,86` | Hard-coded index selection `idx([1,2,5,6,9,10])` for OFF flashes and column deletion for ON — assumes exactly 3 reps with a specific interleaving. If protocol changes, this breaks silently. |
| **Voltage scaling hardcoded** | Multiple files | `v_data = Log.ADC.Volts(2, :)*10` — the ×10 gain factor is repeated in 5+ files rather than centralised. |

### 6.3 Low Priority / Style

| Issue | Detail |
|-------|--------|
| **Scripts vs functions** | `flash4pix_analysis.m`, `quality_analysis.m`, `plot_bar_responses.m` use `close all; clear;` — script-style with workspace pollution. Should be functions. |
| **Debug print statements** | `assess_raw_rec_quality.m:25,27` (`disp(mean_vm)`, `disp(std_vdata)`), `plot_bar_responses.m:51` (`disp(date_str)`) — likely debug leftovers. |
| **Commented-out CSV export** | `extract_metadata_from_folders.m:139–142` — CSV writing code is commented out; metadata extraction has no persistent output format. |
| **Duplicate plotting functions** | `plot_timeseries_polar_bars.m` and `plot_timeseries_polar_bars_pharma.m` share ~90% identical code (pharma version adds drug condition). Could be unified with a flag. |

---

## 7. Cross-Platform Portability

The codebase was primarily developed on **Windows** (G4 arena control) with analysis also run on **macOS**:

- `generate_protocol2.m` and `process_protocol2.m` use Windows `C:\matlabroot\...` paths
- `process_protocol2_pharma.m` has both Windows (active) and Mac (commented) paths
- `plot_bar_responses.m` and `quality_analysis.m` use Mac paths with a different user (`burnettl`)
- No centralised configuration file for paths — each script defines its own `PROJECT_ROOT`

**Recommendation:** Create a shared `config.m` or `get_project_root.m` function that resolves paths based on `computer` or environment variables.

---

## 8. Pipeline Maturity Assessment

| Component | Maturity | Notes |
|-----------|----------|-------|
| Stimulus generation | **Stable** | Well-structured functions with argument validation. Clear separation of patterns vs position functions. |
| Protocol generation | **Stable** | `generate_protocol2.m` is the main entry point with Google Sheets logging. |
| Single-experiment analysis | **Mature** | Two implementations: older `process_protocol2.m` and newer `analyze_single_experiment.m` (cleaner, with `opts` struct and `onCleanup`). |
| Batch analysis | **Active development** | `batch_analyze_1DRF.m` is recent, well-documented, supports median/MAD and mean/SEM statistics. |
| Population analysis | **Mixed** | `combine_bar_results.m` handles both old struct and new table formats. Some population plotting functions (`polar_mean_by_strain.m`, etc.) appear functional. |
| Quality analysis | **Early** | `quality_analysis.m` is a hardcoded script. `assess_raw_rec_quality.m` has debug prints. |
| Pharmacology pipeline | **Parallel fork** | `process_protocol2_pharma.m`, `parse_bar_data_pharma.m`, `plot_*_pharma.m` — duplicated code rather than parameterised. |

---

## 9. Summary of Recommendations

1. **Centralise path configuration** — Replace hardcoded `PROJECT_ROOT` definitions with a single `get_project_root()` function or a config file.
2. **Protect `cd()` calls** — Wrap `load_protocol2_data.m` with `onCleanup` (as `analyze_single_experiment.m` already does) or refactor to use `fullfile()` without `cd()`.
3. **Add error handling to file I/O** — Check `fopen` return values in `save_function_G4_test.m`; capture exceptions in `catch ME` blocks.
4. **Centralise voltage scaling** — Define `VOLTAGE_GAIN = 10` once and reference it, rather than repeating `*10` in every file.
5. **Resolve or remove TODOs** — The `parse_flash_data.m:218` TODO is garbled and should be clarified or deleted.
6. **Unify pharma/standard pipelines** — Refactor `_pharma` variants to accept a drug/condition flag instead of maintaining parallel code.
7. **Convert scripts to functions** — `flash4pix_analysis.m`, `quality_analysis.m`, `plot_bar_responses.m` should be functions to avoid workspace pollution.
8. **Clean up commented-out code** — Rely on git history instead of keeping 100+ lines of dead code in active files.
9. **Add a results export step** — Consider adding CSV/Excel export for final aggregated tables to support downstream analysis outside MATLAB.
10. **Standardise figure output** — All newer code uses `exportgraphics` with vector PDFs, which is good. Ensure older plotting functions follow the same pattern.
