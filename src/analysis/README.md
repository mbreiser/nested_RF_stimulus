# src/analysis

Turns a recorded protocol 2 experiment into direction-selectivity metrics and receptive-field
fits.

## Entry points

| Run this | What it does | Calls |
|---|---|---|
| `process_protocol2()` | Full analysis of one experiment folder | `process_bars_p2`, `process_flash_p2`, `process_bar_flashes_p2` |
| `assess_DS_metrics.m` | Summarises direction-selectivity metrics across cells | — |
| `assess_raw_rec_quality.m` | Recording-quality check before analysis | — |
| `combine_bar_results.m` | Pools bar results across experiments | — |

Run `process_protocol2()` from **inside** the experiment directory, not from the repo root.

## Inputs and outputs

**Reads** (all inside the experiment folder):

- `currentExp.mat` — metadata: frame, age, strain, arena side
- `Log Files/G4_TDMS_Log*.mat` — voltage and frame data
- `Patterns/` and `Functions/` — as written by `generate_protocol2()`

**Writes**, under the project root:

| Path | Contents |
|---|---|
| `results/bar_results/bar_results_*.mat` | DS metrics per cell |
| `results/flash_results/rf_results_*.mat` | RF parameters per cell |
| `figures/bar_stimuli/*.pdf` | Polar plots, heatmaps, timeseries |
| `figures/flash_stimuli/*.pdf` | RF heatmaps, Gaussian fits, contours |

## What each analysis extracts

| Analysis | Detail |
|---|---|
| Bars | 16 directions at 3 speeds (28, 56, 168 deg/s); DSI; preferred direction by vector sum |
| Flashes | 4 px and 6 px grids; 2-D Gaussian RF fit; excitatory and inhibitory components separated |
| Bar flashes | 2 speeds (80 ms, 14 ms); 8 orientations x 11 positions |

## Subdirectories

| Directory | Contents |
|---|---|
| `protocol2/` | Per-analysis scripts for protocol 2 data |
| `analyse_bar_DS/` | Bar direction-selectivity analysis |
| `results_analysis/` | Table building across cells |
| `quality_check/` | Recording and fit quality comparisons |
| `plotting/`, `stats/`, `helper/` | Shared plotting, statistics and utilities |

## Gotchas

- Several scripts in `protocol2/` and `results_analysis/` still contain hardcoded data paths
  near the top; edit them before running on another machine.
