# CLAUDE.md — nested_RF_stimulus

## Project overview

MATLAB codebase for analyzing **Drosophila T4/T5 neuron electrophysiology** using a G4 LED arena system. Two-protocol receptive field (RF) mapping: Protocol 1 (coarse RF) identifies the cell, Protocol 2 (high-res RF + direction selectivity) provides the main analysis data.

**Data location:** `/Users/reiserm/Documents/ttl_1DRF/` (25 experiment folders)
**Default test experiment:** `2025_11_10_10_17`

## What we've done so far

### 1. Codebase audit
- Created `CODEBASE_REVIEW.md` and `ANALYSIS_PIPELINE_REVIEW.md`
- Confirmed only `.mat` files are needed from remote (no `.tdms`)

### 2. Single experiment validation
- Created `scripts/validate_and_analyze_single.m` for one-shot testing
- Fixed `circ_vmpar` error by adding CircStat toolbox to path

### 3. Enhanced visualization (`analyze_single_experiment_mr.m`)
Created `src/analysis/protocol2/analyze_single_experiment_mr.m` — a copy-and-extend variant of `analyze_single_experiment.m` with `opts.visualize_everything` enhancements:

| Enhancement | Description | Status |
|---|---|---|
| A | Voltage histogram figure (`plot_voltage_histograms.m`) | Working |
| B | Green stimulus timing lines on bar flash figures | Working on Figs 2 & 4; **broken on Fig 3** |
| C | PD row highlight on 8x11 heatmap | Working |
| D | PD (red) + orthogonal (black) direction lines on polar plot | Working |
| E | Relabel ND/PD to Leading/Trailing on Fig 3 | Working |
| F | Time scale bars on all figures | Working |

**New files created (not modifying originals):**
- `src/analysis/protocol2/analyze_single_experiment_mr.m`
- `src/analysis/plotting/plot_voltage_histograms.m`
- `src/analysis/plotting/add_stim_timing_lines.m`

### Known bug: Fig 3 green timing lines
The green stimulus timing lines disappear from Figure 3 (1x11 PD-ND bar flash plot). Root cause: MATLAB tiledlayout `Position` queries trigger layout recalculations that reset the `Children` stack reordering used for z-ordering. We moved Enhancement F (scale bars) before Enhancement B (timing lines) to mitigate, but the issue persists on Fig 3 specifically. The lines are programmatically present (verified via handle inspection) but not visible in the rendered output or saved PDFs. Needs further investigation — may require an alternative approach to z-ordering (e.g., drawing lines first, then data on top, or using `uigridlayout`).

## Key technical details

- **Sampling rate:** 10 kHz (1 sample = 0.1 ms)
- **Voltage scaling:** `Log.ADC.Volts(2,:) * 10` converts to mV
- **Frame signal:** `Log.ADC.Volts(1,:)` encodes stimulus identity; 0 = grey screen
- **Bar sweeps:** 16 directions x 3 speeds (28, 56, 168 dps) x 3 reps — only 28 dps plotted
- **Bar flashes:** 8 orientations x 11 positions — 80ms (slow) and 14ms (fast) — only 80ms plotted
- **Flash timing:** onset = sample 5001, offset = sample 5801, duration = 801 samples (80.1ms)
- **CircStat toolbox:** `/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a/`

### Quick run commands
```matlab
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

exp_folder = '/Users/reiserm/Documents/ttl_1DRF/2025_11_10_10_17';
opts = struct();
opts.save_figs = true;
opts.save_dir = fullfile(exp_folder, 'analysis_output_mr');
opts.visualize_everything = true;
analyze_single_experiment_mr(exp_folder, opts);
```

## Pipeline architecture: single fly to population

```
Single experiment folder
    ↓
analyze_single_experiment.m (or _mr variant)
    ↓  parse_bar_data → compute_bar_sweep_responses → find_pd_from_lut
    ↓  parse_bar_flash_data → extract per-cell traces
    ↓
batch_analyze_1DRF.m  (loops all experiments)
    ↓  process_single_cell() for each folder
    ↓  classify: ON/OFF (frame > 129) × control/TTL (strain name)
    ↓  extract_flash_traces() → 11×N baseline-subtracted matrices
    ↓  compute_peak_metrics() → peak_pos (M2), peak_amplitudes, bump_width
    ↓  compute_m5_centroid() → centroid_m5, centroid_m5_rounded
    ↓  reindex_to_peak() → M2-aligned 11×N traces (M2 peak at row 6)
    ↓  reindex_to_peak() → M5-aligned 11×N traces (round(M5) at row 6)
    ↓  extract_temporal_metrics() → rise/decay timing (M2-based + M5-based)
    ↓  voltage metrics: rest/stim median/mean
    ↓
Population plots (per ON and OFF group):
    - plot_polar_population.m (PD-aligned tuning curves, mean±SEM)
    - plot_flash_1x11_population.m (unaligned, M2-aligned, M5-aligned)
    - plot_baseline_voltage_comparison.m (resting vs stim Vm by group)
    - plot_rf_width_comparison.m (FWHM bump width by group)
    - plot_dsi_comparison.m (DSI by group)
    - plot_ortho_width_comparison.m (PD vs orthogonal RF width)
    - plot_timing_by_position.m (temporal metrics vs position, M2 + M5)
    - plot_amplitude_by_position.m (per-position depol/hyperpol, rank-sum)
    - plot_response_histograms_by_position.m (per-position response distributions)
    - plot_tuning_t4_vs_t5.m (T4 vs T5 FWHM and DSI comparison)
```

**Output:** Population figures + preview PNGs saved to `<data_root>/population_results/` and `<data_root>/figure_previews/`

### Quick run: generate all preview PNGs
```matlab
run('scripts/generate_figure_previews.m')       % main pipeline + polar/flash/timing PNGs
run('scripts/run_position_distributions.m')      % per-position amplitude + histogram PNGs
run('scripts/run_t4_vs_t5_comparison.m')         % T4 vs T5 polar + tuning PNGs
```

## Baseline definitions — track carefully

Different parts of the pipeline use different baseline windows. This needs validation to determine whether the differences matter for final analysis.

| Context | Baseline window | Duration | Code location |
|---|---|---|---|
| Bar sweep responses | samples 1000–9000 | 0.8 s | `compute_bar_sweep_responses.m` |
| Bar flash (single cell plots) | samples 1–5000 | 0.5 s | `plot_bar_flash_1x11.m` |
| Bar flash (batch extraction) | samples 1–5000 | 0.5 s | `batch_analyze_1DRF.m:extract_flash_traces()` |
| Bar flash (population plots) | pre-subtracted | — | `plot_flash_1x11_population.m` |

The bar flash baseline (1:5000) covers the full pre-stimulus period (stim onset at sample 5001). The bar sweep baseline (1000:9000) skips the first 1000 samples and uses a longer window. Both should give similar results in practice but the sweep baseline intentionally avoids edge artifacts from the trace extraction.

**Open question:** Does the sweep baseline choice (excluding first 1000 samples) meaningfully affect PD computation or peak amplitude? This should be checked by comparing results with a 1:9000 baseline.

## Upstream sync (2025-02-17)

Synced local repo with Laura's latest `origin/main` (8 commits). Our work lives on the `mr-enhancements` branch, rebased cleanly on top.

**What Laura changed:** Removed pharma-related code (5 plotting functions, `parse_bar_data_pharma.m`, `process_bars_p2_pharma.m`). Added new helpers (`addOrthoMetrics.m`, `shiftMaxColumnTo5.m`, `runGroupedStats.m`, `plotGroupedBox.m`, `plotPolarByGroup.m`). Added stimulus schematic figures and experiment log spreadsheet.

**API break handled:** `parse_bar_flash_data` now requires a 3rd argument `prop_int` (proportion of inter-flash interval to keep). Old hardcoded `gap_between_flashes = 5000` equals `prop_int = 0.5`. Laura's newer `process_bar_flashes_p2.m` uses `prop_int = 0.75`. Our `_mr` variant uses `0.5` to match old behavior and keep timing calculations consistent (stim onset at sample 5001 with gap=5000).

**Note:** `analyze_single_experiment.m` still calls `parse_bar_flash_data` with only 2 args — it'll error on the updated code. `batch_analyze_1DRF.m` has been fixed to pass `prop_int = 0.5`. This is Laura's incomplete migration; flag for her regarding the original `analyze_single_experiment.m`.

### 4. RF centroid algorithm overhaul

Rewrote `compute_rf_centroid` in `analyze_single_experiment_mr.m` to use a bump-first approach:
- **Amplitude**: 99.5th percentile in response window (samples 5001–6551), baseline-subtracted, negatives thresholded to 0
- **Bump finding**: FWHM — contiguous region around peak where amplitude ≥ peak/2
- **Centroid**: response-weighted center-of-mass over bump positions only
- Returns `rf_metrics` struct with peak_pos, peak_val, bump_range, bump_width, centroid_peak_delta

Also updated `annotate_centroid` to show bump tiles (blue borders) and centroid tile (red border) on Fig 3.

**New files:**
- `scripts/compare_centroid_methods.m` — old vs new centroid comparison (25 cells)
- `scripts/compare_five_centroid_methods.m` — 5-method comparison (see below)
- `scripts/regenerate_all_centroid_figures.m` — batch regeneration driver

**Modified:**
- `src/analysis/protocol2/analyze_single_experiment_mr.m` — rewritten `compute_rf_centroid`, updated caller & output struct, rewritten `annotate_centroid`
- `scripts/run_on_cells_mr.m` — updated summary format (Peak@, Bump, Delta columns)
- `scripts/run_off_cells_mr.m` — same

### 5. Five-method centroid comparison

Compared 5 RF center-finding methods on all 25 cells (`scripts/compare_five_centroid_methods.m`):

| Method | Description | Output |
|---|---|---|
| M1 | Assumed center (always pos 6) | 6.00 (constant) |
| M2 | Robust peak (position of max amplitude) | integer |
| M3 | Flanked peak (M2 + validation: both neighbors must rank in top 4) | integer or "no bump" |
| M4 | Full depolarizing centroid (center-of-mass over all A>0 positions) | fractional |
| M5 | FWHM bump centroid (center-of-mass restricted to half-max region) | fractional |

**Key results:**
- 13/25 cells flagged (MaxΔ > 1.0 position between M2–M5, or M3 = no bump)
- 6/25 cells fail M3 flanked-peak test
- M4 vs M5 mean difference small (0.07–0.27 positions) but M4 pulls toward center
- Bump width by group: ON ctrl 5.0±0.7, ON TTL 7.2±2.3, OFF ctrl 5.3±1.3, OFF TTL 7.0±1.3
- TTL cells have broader bumps — potential phenotype

### RF alignment decision (2025-02-19)

**Initial decision: Use M2 (peak position / argmax) for population alignment.** This matches Gruntman et al. (eLife 2019) and avoids the differential asymmetry bias between control (narrow bumps) and TTL (broad bumps) that M5 introduces. Integer resolution only, but unbiased.

**Update (2025-02-20):** Both M2 and M5 alignment are now run in parallel through the full population pipeline (see §7 below). This allows direct visual comparison to inform the final alignment choice. 13/25 cells have different M2 vs round(M5) positions.

### 6. Gruntman alignment + population pipeline update (2025-02-19)

Aligned the batch pipeline with Gruntman et al. methods and added new per-cell metrics.

**Key decisions (user-confirmed):**
- RF alignment: M2 (peak position) — matches Gruntman
- Population stats: mean ± SEM — matches Gruntman
- Peak detection: keep 99.5th percentile (deviation from Gruntman's max())

**Changes to `batch_analyze_1DRF.m`:**
- Default stat_method changed from `'median_mad'` to `'mean_sem'`
- Added per-cell RF metrics: `peak_pos`, `peak_amplitudes`, `bump_width`
- Added per-cell voltage metrics: `rest_voltage_{median,mean}`, `stim_voltage_{median,mean}`
- Added peak-aligned traces: `pd_flash_peak_aligned`, `ortho_flash_peak_aligned`
  - Uses `reindex_to_peak()`: shifts 11×N traces so peak_pos maps to row 6 (center)
  - Off-edge rows become NaN (handled by omitnan in population stats)
- Added temporal metrics: `temporal_metrics` struct via `extract_temporal_metrics()`
- Added peak-aligned population plot calls with `require_both_n2 = true`
- New local functions: `compute_peak_metrics`, `reindex_to_peak`, `compute_bump_width`, `save_aligned_figures`

**Changes to `plot_flash_1x11_population.m`:**
- Added `opts.require_both_n2` option (default: false)
- When true, positions are only plotted if BOTH ctrl and TTL have n≥2

**New files:**
- `src/analysis/helper/extract_temporal_metrics.m` — Gruntman-style position-by-position timing (rise start 10%, rise time 10→50%, decay time 80→20%)
- `src/analysis/plotting/plot_baseline_voltage_comparison.m` — resting vs stimulus voltage by group
- `src/analysis/plotting/plot_rf_width_comparison.m` — RF FWHM width by group with Wilcoxon tests
- `src/analysis/plotting/plot_timing_by_position.m` — temporal metrics vs position (mean±SEM)
- `scripts/generate_example_cell_panels.m` — scaffold for example cell sweep/polar PNGs (user fills in folder names)
- `scripts/generate_figure_previews.m` — master PNG generation script for all population figures

### 7. Dual M2/M5 alignment comparison (2025-02-20)

Extended the population pipeline to run **both M2 (argmax) and M5 (FWHM bump centroid)** alignment in parallel, producing side-by-side labeled summary plots for visual comparison.

**M5 fractional → integer:** `round(centroid_m5)`, clamped to [1,11]. 13/25 cells have different M2 vs round(M5) positions.

**New per-cell fields in `batch_analyze_1DRF.m`:**
- `centroid_m5` — fractional M5 centroid (e.g. 5.55)
- `centroid_m5_rounded` — integer for alignment (e.g. 6)
- `m5_bump_range` — [left, right] of FWHM bump
- `pd_flash_m5_aligned`, `ortho_flash_m5_aligned` — M5-aligned 11×N traces
- `temporal_metrics_m5` — Gruntman-style timing computed from M5 alignment

**New local function:** `compute_m5_centroid(A)` — FWHM walk-from-peak algorithm (same as `compute_rf_centroid` in `analyze_single_experiment_mr.m`), returns struct with `.centroid`, `.centroid_int`, `.bump_range`, `.bump_width`.

**Changes to `plot_flash_1x11_population.m`:**
- Added `opts.resp_end_sample` — draws thin dark vertical line at response window end (sample 6551) on each tile
- Added `opts.show_ordinal_ranks` — annotates ordinal amplitude rank (#1=strongest, #2=next, etc.) in upper-right corner of each tile, separately for ctrl (black) and ttl (red)
- Restructured main loop into pre-compute pass (stats + ranks) then plot pass

**Changes to `plot_timing_by_position.m`:**
- Added `opts.temporal_field` (default: `'temporal_metrics'`) — callers pass `'temporal_metrics_m5'` for M5

**Changes to `generate_figure_previews.m`:**
- Completely rewritten for dual M2/M5 output
- Renamed `*_peak_aligned.png` → `*_M2_aligned.png`
- Added M5-aligned PNGs and M5 timing PNGs

**Preview PNGs generated (22 total in `<data_root>/figure_previews/`):**

| Category | PNGs |
|---|---|
| Polar tuning (2) | `pop_polar_{on,off}_mean_sem.png` |
| Unaligned flash (4) | `pop_{pd,ortho}_flash_{on,off}_mean_sem.png` |
| M2-aligned flash (4) | `pop_{pd,ortho}_flash_{on,off}_M2_aligned.png` |
| M5-aligned flash (4) | `pop_{pd,ortho}_flash_{on,off}_M5_aligned.png` |
| Summary comparisons (4) | `pop_baseline_voltage.png`, `pop_rf_width.png`, `pop_dsi_comparison.png`, `pop_ortho_width_comparison.png` |
| Timing M2+M5 (4) | `pop_timing_{on,off}_{M2,M5}.png` |

**Aligned flash plot annotations (on all M2 and M5 aligned plots):**
- Green vertical lines at samples 5001 (stim onset) and 5801 (stim offset)
- Thin dark vertical line at sample 6551 (response window end / peak detection boundary)
- Ordinal amplitude rank labels (#1–#11) per group per position
- Per-position n-count labels
- Custom center labels: "0 (M2 peak)" or "0 (M5 centroid)"
- `require_both_n2 = true` (positions only plotted if both groups have n≥2)

### 8. Additional population metrics (2025-02-20)

Added DSI, direction tuning width, orthogonal RF metrics to the batch pipeline and individual cell figures.

**New per-cell fields in `batch_analyze_1DRF.m`:**
- `dsi` — direction selectivity index
- `tuning_width_deg` — direction tuning width (degrees, from von Mises fit)
- `ortho_bump_width` — RF width along orthogonal axis
- `pd_ortho_aspect_ratio` — PD bump width / orthogonal bump width

**New plotting functions:**
- `src/analysis/plotting/plot_dsi_comparison.m` — DSI by group with rank-sum stats
- `src/analysis/plotting/plot_ortho_width_comparison.m` — PD vs orthogonal RF width comparison

**Individual cell enhancements:**
- Polar plot annotations now show DSI, tuning width, ortho bump width, and PD/ortho aspect ratio
- All 25 individual cell figures regenerated with new annotations

### 9. Time-to-90% metric and plot refinements (2025-02-23)

Added a simple, intuitive temporal metric and improved aligned flash plot annotations.

**New temporal metric: `time_to_90`**
- Time from flash onset (sample 5001) to first sample reaching 90% of peak response
- Peak = strict `max()` of mean trace (averaged across 3 reps) — appropriate for timing on smoothed data
- Added to `extract_temporal_metrics.m` as both raw (`time_to_90`) and delta (`time_to_90_delta`) fields

**Timing plots switched from delta to raw values:**
- All 4 panels in `plot_timing_by_position.m` now show raw (absolute) metrics instead of delta-from-peak
- **Reason:** Delta normalization subtracts each cell's peak-position value independently, forcing both ctrl and TTL to ~0 at center and washing out group timing differences
- Raw values preserve the absolute timing comparison between groups (e.g., ctrl ~80ms vs TTL ~90ms at center)
- Fields plotted: `rise_start`, `rise_time`, `decay_time`, `time_to_90`

**Stimulus timing lines on aligned flash population plots:**
- Added green vertical lines at stim onset (sample 5001) and stim offset (sample 5801) to `plot_flash_1x11_population.m`
- New options: `opts.stim_onset_sample`, `opts.stim_offset_sample`
- Three-line scheme: two green (stim boundaries) + one dark gray (response window end at 6551)
- Applied to all M2 and M5 aligned flash plots via both `batch_analyze_1DRF.m` and `generate_figure_previews.m`

**Files modified:**
- `src/analysis/helper/extract_temporal_metrics.m` — added `time_to_90` metric
- `src/analysis/plotting/plot_timing_by_position.m` — expanded to 4x1 layout, switched to raw values
- `src/analysis/plotting/plot_flash_1x11_population.m` — added green stim timing lines
- `src/analysis/protocol2/batch_analyze_1DRF.m` — passes stim timing samples to aligned flash opts
- `scripts/generate_figure_previews.m` — passes stim timing samples to aligned flash opts

**Note on batch_results.mat save gate:** The results .mat file save in `batch_analyze_1DRF.m` is controlled by `opts.save_figs` (line 142), not a separate `save_results` flag. Must set `opts.save_figs = true` to trigger the save.

### Gruntman discrepancy summary

| Aspect | Gruntman | Ours | Decision |
|---|---|---|---|
| Sampling rate | 20 kHz | 10 kHz | Keep — Nyquist-limited |
| Flash baseline | 900ms pre-stim | 500ms (1:5000) | Keep — full pre-stim period |
| Peak detection | max() | 99.5th percentile | Keep ours |
| RF alignment | Peak position (argmax) | M2 (argmax) | Matched |
| Pop stats | mean + SEM | mean + SEM | Matched |
| Timing metrics | rise/decay per position (delta) | Implemented (raw values) | Diverged — raw preserves group differences |
| Trial exclusion | ±10/15 mV thresholds | Not implemented | Deferred (only 3 reps) |

### 10. Orthogonal direction fix — right-hand rule (2025-02-23)

Fixed a critical alignment bug: 12 of 25 cells (48%) had ortho bar flash positions flipped relative to the right-hand rule convention.

**Problem:** `find_pd_from_lut.m` computed a single `pos_order` from PD function parity (odd=forward, even=reverse) and applied it to **both** PD and orthogonal traces. The orthogonal bar's forward sweep direction has no consistent relationship to the PD function parity, so nearly half the cells had inverted ortho positions.

**Solution — right-hand rule:**
1. Look up ortho bar's forward sweep direction from LUT: `Tbl.direction` for `function==3`
2. Compute right-hand target: `rh_target = mod(pd_direction - 90, 360)` (90° CW from PD)
3. Angular distance between `ortho_dir_fwd` and `rh_target`: if < 90° → `ortho_pos_order = 1:11`, else → `11:-1:1`

**Diagnostic:** `scripts/diagnose_ortho_direction.m` — verified all 25 cells. 13 match, 12 need flipping. Generates polar diagnostic figure (`diag_ortho_direction.png`).

**Files modified:**
- `src/analysis/helper/find_pd_from_lut.m` — added `ortho_pos_order` via right-hand rule, new output field `pd_info.ortho_pos_order`
- `src/analysis/protocol2/batch_analyze_1DRF.m` — ortho trace extraction (`ortho_flash_bl`) and ortho peak metrics (`compute_peak_metrics_for_col`) now use `pd_info.ortho_pos_order` instead of `pd_info.pos_order`

**New files:**
- `scripts/diagnose_ortho_direction.m` — diagnostic script for verifying ortho direction geometry

**Note:** This change requires re-running the batch pipeline to regenerate `batch_results.mat` with corrected ortho traces. PD traces are unaffected.

### 11. Per-position amplitude distributions (2025-02-23)

Added two new population analysis functions for per-position amplitude metrics with statistical tests.

**New files:**
- `src/analysis/plotting/plot_amplitude_by_position.m` — 2×1 figure: depolarization (99.5th pctile) and hyperpolarization (0.5th pctile) at each aligned position, jittered dots + mean±SEM, Wilcoxon rank-sum asterisks
  - Rejection criterion: positions where `|mean(trace(5001:end))| < threshold` are excluded
  - Expanded hyperpolarization window: `[5001, end]` (vs `[5001, 6551]` for depolarization)
- `src/analysis/plotting/plot_response_histograms_by_position.m` — 2×6 grid of histogram panels, each showing per-cell mean response distribution (stim onset to end of trace) for ctrl vs TTL, with rank-sum p-values
- `scripts/run_position_distributions.m` — driver script generating 8 PNGs

**PNGs generated (8):**
- `pop_amp_dist_{on,off}_{M2,M5}.png` (4 — per-position peak metrics)
- `pop_resp_hist_{on,off}_{M2,M5}.png` (4 — response distributions)

### 12. T4 vs T5 tuning comparison (2025-02-23)

Added T4 (ON) vs T5 (OFF) tuning comparison using blue=T4, orange=T5 color scheme.

**Modified:** `src/analysis/plotting/plot_polar_population.m`
- Added optional group label/color overrides: `opts.group1_label`, `opts.group2_label`, `opts.group1_line_color`, `opts.group1_fill_color`, `opts.group2_line_color`, `opts.group2_fill_color`
- Default behavior unchanged (ctrl=black, TTL=red)

**New files:**
- `src/analysis/plotting/plot_tuning_t4_vs_t5.m` — 1×2 figure: direction tuning FWHM and DSI, comparing T4 vs T5 within each treatment (ctrl and TTL), with Wilcoxon rank-sum brackets
- `scripts/run_t4_vs_t5_comparison.m` — driver script generating 3 PNGs

**PNGs generated (3):**
- `pop_polar_ctrl_t4_vs_t5.png` — Polar overlay: T4-ctrl vs T5-ctrl (blue vs orange)
- `pop_polar_ttl_t4_vs_t5.png` — Polar overlay: T4-TTL vs T5-TTL (blue vs orange)
- `pop_tuning_t4_vs_t5.png` — FWHM + DSI box plots: T4 vs T5

### 13. Pre-bar-flash dataset — directional tuning batch analysis (2025-02-26)

Analyzed a second dataset from an **earlier protocol** (summer 2025) at `/Users/reiserm/Documents/ttl_1DRF/pre-bar-flash/`. These experiments have bar sweeps at 2 speeds (28 + 56 dps) but **no bar flashes**, no 168 dps speed, and no RF mapping data. Pipeline extracts directional tuning and computes per-cell DS metrics.

**Dataset:** 23 experiments in `{control,ttl}/{ON,OFF}/` directory tree:

| Group | n | Notes |
|-------|---|-------|
| ON control | 7 | Includes 1 post-Oct-15 (`2025_10_20_11_53`) |
| OFF control | 6 | Includes 1 post-Oct-15 (`2025_10_22_13_10`) |
| ON TTL | 8 | All pre-Oct-15 |
| OFF TTL | 2 | All pre-Oct-15 |

**Key design decisions:**
- **No modifications to existing code** — all new files, reuses existing functions unchanged
- **Classification by directory path**, not metadata — post-Oct-15 OFF control experiment has wrong `metadata.Strain = '42F06_T4T5_ttl'`; directory position is the ground truth
- **Duration-based parser** — original `parse_bar_data.m` hardcodes gap indices for 3-speed + bar-flash protocol. New parser classifies bar segments by physical properties (frame std > 5 for moving, duration > 1.5s for slow/28dps) to automatically select only 28dps bars regardless of protocol version
- **Dark-bar polarity correction** — see §13a below

**New files:**

| File | Description |
|------|-------------|
| `scripts/diagnose_pre_bar_flash.m` | Diagnostic — validated gap structure, segment classification, LUT directions on sample experiments |
| `scripts/diagnose_off_direction_swap.m` | Diagnostic — compared LUT direction mappings between ON and OFF cells |
| `scripts/diagnose_off_direction_swap_v2.m` | Diagnostic — analyzed frame signal sweep direction per slow bar segment |
| `src/analysis/protocol2/pipeline/parse_bar_data_pre_bf.m` | 2-speed bar data parser (duration-based segment classification) |
| `src/analysis/protocol2/batch_analyze_pre_bar_flash.m` | Batch pipeline — bar sweeps only, no bar flashes or RF mapping |
| `src/analysis/helper/correct_off_polarity_swap.m` | Dark-bar polarity correction for pre-Oct-15 OFF cells (§13a) |
| `scripts/run_pre_bar_flash_analysis.m` | Driver — runs batch + generates population plots + metrics table |
| `scripts/plot_pre_bf_sweeps_individual.m` | Per-cell polar ring-of-traces figures with DSI annotation |

**Reused functions (unchanged):** `load_protocol2_data`, `compute_bar_sweep_responses`, `find_PD_and_order_idx`, `compute_bar_response_metrics`, `verify_lut_directions`, `plot_polar_population`, `plot_dsi_comparison`, `plot_aspect_ratio_comparison`

**Parser algorithm (`parse_bar_data_pre_bf.m`):**
1. Find cycle boundaries from grey-screen gaps ≥ 8s (exactly 4 gaps → 3 cycles)
2. Per cycle: find all frame transitions (`|diff(f)| > 9`), build segments between transitions
3. Classify segments: moving (frame std > 5) AND slow (duration > 15000 samples = 1.5s @ 10kHz)
4. Extract 16 slow-bar traces per cycle with 9000-sample pre/post padding
5. Average across 3 cycles → 16×4 cell array matching `parse_bar_data` output contract

**Per-cell results struct fields:** `folder`, `date_str`, `strain`, `frame`, `is_on`, `is_ttl`, `group`, `max_v_aligned` (16×2), `pd_direction`, `dsi_vector`, `dsi_pdnd`, `dir_tuning_fwhm`, `dir_tuning_cv`, `dir_tuning_kappa`, `tuning_aspect_ratio`, `sym_ratio`, `polarity_corrected`

**Key results (after polarity correction):**

| Group | n | DSI_vector | DSI_pdnd | FWHM | Aspect Ratio |
|-------|---|-----------|----------|------|-------------|
| ON control | 7 | 0.241±0.037 | 0.392±0.096 | 212±101° | 2.78±0.44 |
| ON TTL | 8 | 0.200±0.048 | 0.311±0.105 | 307±56° | 1.68±0.30 |
| OFF control | 6 | 0.230±0.028 | 0.328±0.048 | 304±53° | 2.67±0.50 |
| OFF TTL | 2 | 0.227±0.063 | 0.261±0.096 | 338±0° | 2.07±0.28 |

ON TTL aspect ratio (1.68) substantially lower than ON control (2.78), consistent with broader tuning under TTL. OFF TTL has only n=2 — insufficient for statistics.

**Output PNGs (4) + metrics table in `<data_root>/figure_previews/`:**
- `pre_bf_polar_on.png`, `pre_bf_polar_off.png` — PD-aligned polar tuning (mean±SEM)
- `pre_bf_dsi_comparison.png` — DSI + FWHM by group
- `pre_bf_aspect_ratio.png` — Aspect ratio: T4 vs T5 × treatment
- `pre_bf_metrics_table.txt` — Per-cell metrics with summary statistics

**Quick run:**
```matlab
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');
run('scripts/run_pre_bar_flash_analysis.m');
```

**Saved results:** `/Users/reiserm/Documents/ttl_1DRF/pre-bar-flash/population_results/batch_results_pre_bf.mat`

### 13a. Dark-bar polarity correction for pre-Oct-15 OFF cells (2025-02-26)

**Problem:** All pre-Oct-15 OFF cells had 2 pairs of directions visibly swapped on polar tuning curves, while ON cells looked correct.

**Root cause:** Pattern files `0009` and `0010` (repo: `results/patterns/protocol2/full_field_bars4/`) use `bar_0` (dark bar, pixel values 0–3 on background=4) instead of `bar_15` (bright bar, pixel values 5–15 on background=4). Commit `82badc7c` ("update 15 and 16 patterns — correct direction", Oct 15 2025) fixed patterns 0015/0016 from `bar_0` → `bar_15`, but patterns 0009–0014 remain with `bar_0` in the repository. The arena SD card patterns were updated around the same time.

**What the polarity inversion does:**
- The bar sweep direction (motion angle) is **unchanged** — bright and dark bars traverse the same path
- However, for OFF cells (which respond to dark features), the polarity inversion causes the func=3 and func=4 response amplitudes to appear **swapped** for patterns 9/10
- ON cells are unaffected because their excitatory response direction is the same for bright vs dark bars
- Experiment patterns 9 and 10 map to 4 of 16 slow bar directions: 22.5°, 45°, 202.5°, 225°

**Empirical verification (all 23 cells):**
- 15/15 ON cells: original data correct (swap worsens smoothness metric)
- 7/7 pre-Oct-15 OFF cells: swap dramatically improves smoothness
- 1/1 post-Oct-15 OFF cell (`2025_10_22_13_10`): original data correct (patterns already fixed on arena)

**Correction implementation (`correct_off_polarity_swap.m`):**
1. Check if experiment is OFF cell AND date < Oct 15, 2025
2. If so, swap bar_data rows between func=3 and func=4 for patterns 9 and 10
3. Applied automatically in both `batch_analyze_pre_bar_flash.m` and `plot_pre_bf_sweeps_individual.m`
4. Result stored in `results.polarity_corrected` (true/false) for traceability

**Impact on results:**
- OFF control DSI_vector: 0.077±0.057 → **0.230±0.028** (now comparable to ON control 0.241)
- OFF control aspect ratio: 1.43±0.51 → **2.67±0.50** (now comparable to ON control 2.78)
- OFF TTL DSI_vector: 0.081±0.004 → **0.227±0.063**
- ON cells: **unchanged** (correction not applied)

**Diagnostic scripts used:**
- `scripts/diagnose_off_direction_swap.m` — LUT mapping comparison between ON and OFF cells
- `scripts/diagnose_off_direction_swap_v2.m` — frame signal sweep direction analysis per segment

### 14. DS alignment exploration — 4 PD extraction methods (2025-02-27)

Compared 4 PD extraction/alignment methods on the combined dataset (48 cells: 23 early + 25 late) to assess whether finer alignment precision or symmetry-based PD improves population tuning curve shape.

**4 methods:**

| Method | Short name | PD computation | Output grid |
|--------|-----------|----------------|-------------|
| 1 | VecNN | Vector sum → snap to nearest 22.5° grid (current) | 16-point (22.5°) |
| 2 | SymNN | Symmetry axis from 16 candidates → snap to 22.5° | 16-point (22.5°) |
| 3 | VecInterp | Vector sum (continuous) → interpolated 2.5° grid | 144-point (2.5°) |
| 4 | SymInterp | Brute-force symmetry axis at 2.5° → interpolated grid | 144-point (2.5°) |

**New files:**

| File | Description |
|------|-------------|
| `src/analysis/helper/compute_pd_four_methods.m` | Core function: 4 PD methods for one cell's tuning curve. Returns PD angle, aligned tuning curve, FWHM, symmetry score for each method |
| `scripts/ds_alignment_exploration.m` | Driver: loads both batch results, runs 4 methods on all 48 cells, generates individual polar plots, disagreement table, population polars, FWHM/symmetry comparison |

**Modified:** `src/analysis/plotting/plot_polar_population.m` — relaxed size check from `isequal(size(d), [16 2])` to `size(d, 2) == 2 && size(d, 1) >= 2` to accept 144×2 interpolated curves.

**Key results:**
- 6/48 cells flagged with MaxΔ > 22.5° between methods — all due to SymInterp placing PD at 0°/180° (orthogonal axis, not PD axis) in broadly tuned cells
- VecNN always returns 90° because input is already VecNN-aligned `max_v_aligned` data
- SymNN differs from VecNN in 7/48 cells (shifts to 67.5° or 112.5°)
- VecInterp shows fine-grained offsets (typically ±5–10° from 90°)
- Mean MaxΔ by group: ON ctrl 6.3°, ON TTL 22.6°, OFF ctrl 17.3°, OFF TTL 35.1°

**Output PNGs (16 files in `<data_root>/figure_previews/`):**
- `ds_explore_individual_*.png` (6) — per-cell polar with 4 colored PD arrows
- `ds_explore_polar_{VecNN,SymNN,VecInterp,SymInterp}_{on,off}.png` (8) — population polar tuning
- `ds_explore_fwhm_comparison.png` — FWHM by method × group
- `ds_explore_symmetry_comparison.png` — symmetry score by method × group
- `ds_explore_disagreement.txt` — per-cell PD angle table with disagreement flags

**Quick run:**
```matlab
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');
run('scripts/ds_alignment_exploration.m');
```

### 14a. Combined batch comparison (2025-02-27)

Generated combined population plots merging early (23 cells) and late (25 cells) batches.

**New file:** `scripts/run_combined_batch_comparison.m` — loads both batch results, harmonizes into common struct, generates plots for 3 groupings: combined (n=48), early only, late only.

**Output PNGs (12 in `<data_root>/figure_previews/`):**
- `comb_polar_{on,off}_{combined,early,late}.png` (6) — polar tuning
- `comb_dsi_{combined,early,late}.png` (3) — DSI comparison
- `comb_aspect_ratio_{combined,early,late}.png` (3) — aspect ratio

## Figure Layout Composer workflow

We use a visual playground (`tools/figure-layout-composer.html`) to mock up multi-panel figure layouts before generating polished vector PDFs. The workflow is:

1. Generate rough **preview PNGs** of each candidate panel from MATLAB analysis
2. Open the Figure Layout Composer in a browser, select the target journal
3. Drag PNGs onto the canvas, arrange panels, set sizes/positions
4. For each panel, set its **type** (Plot or Image) and write a **description** of what it shows
5. Add annotations (panel-specific or global instructions)
6. Click **Copy Spec** → paste into Claude
7. Claude writes MATLAB/Python code that produces the final vector PDF

**PNGs are throwaway previews only.** The final output is always code → editable vector PDF.

### Preview PNG naming convention

When generating preview PNGs for the layout tool, use this naming scheme so Claude can unambiguously identify which panel is which:

```
fig<N>_<label>_<short_description>.png
```

Examples:
- `fig1_a_polar_tuning_ON_ctrl.png`
- `fig1_b_bar_flash_PD_ND_ON_ctrl.png`
- `fig2_a_population_polar_median_mad.png`
- `fig3_c_voltage_histogram_OFF_ttl.png`

Rules:
- `fig<N>` = which composite figure this belongs to (fig1, fig2, etc.)
- `<label>` = the panel label (a, b, c...) matching the layout composer
- `<short_description>` = plot type + key identifiers (cell group, condition, metric)
- Use underscores, no spaces, all lowercase
- Store previews in `<experiment_folder>/figure_previews/` (single-cell) or `<data_root>/figure_previews/` (population)

### Journal presets available

Nature / Nat Neuro / Nat Methods, Cell / Neuron / Current Biology, Science (AAAS), eLife, PNAS, J Neuroscience, PLOS Biology. Each preset auto-configures canvas size, label formatting, font sizes, and constraint warnings.

## Next steps

1. **Re-run batch pipeline** — Required after ortho direction fix (§10). Run `generate_figure_previews.m` to regenerate `batch_results.mat` with corrected ortho traces, then re-run `run_position_distributions.m` and `run_t4_vs_t5_comparison.m` on the updated results
2. **Review new PNGs** — Inspect the 11 new PNGs (amplitude distributions, response histograms, T4 vs T5 comparisons) plus the ortho diagnostic figure
3. **Review M2 vs M5 alignment** — Visual inspection to decide which alignment method to use going forward
4. **Select example cells** — Fill in folder names in `scripts/generate_example_cell_panels.m` for the 4 representative cells (one per group)
5. **Compose figures** — Use the figure layout composer with the generated PNGs to design the final multi-panel figures
6. **Baseline sensitivity check** — Compare single-cell PD and peak amplitude using different baseline windows to determine if the discrepancy matters
7. **Fix Fig 3 green lines** — Investigate alternative z-ordering approaches
