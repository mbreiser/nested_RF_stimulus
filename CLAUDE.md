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
    ↓
Population plots (per ON and OFF group):
    - plot_polar_population.m (PD-aligned tuning curves)
    - plot_flash_1x11_population.m (bar flash with median±MAD or mean±SEM)
```

**Output:** 6 population figures saved to `<data_root>/population_results/`

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

**Note:** `analyze_single_experiment.m` and `batch_analyze_1DRF.m` still call `parse_bar_flash_data` with only 2 args — they'll error on the updated code. This is Laura's incomplete migration; flag for her or fix locally when running batch analysis.

## Next steps

1. **Validate population pipeline** — Run `batch_analyze_1DRF.m` on the full dataset, verify cell classification (ON/OFF, control/TTL), check that all 25 experiments process without errors. Will need to fix the missing `prop_int` arg in `batch_analyze_1DRF.m` first.
2. **Baseline sensitivity check** — Compare single-cell PD and peak amplitude using different baseline windows to determine if the discrepancy matters
3. **Fix Fig 3 green lines** — Investigate alternative z-ordering approaches
4. **Figure 4 updates** — Orthogonal axis figure may need similar label/annotation updates as Fig 3 (per user, to be decided later)
