# Analysis Pipeline: Step-by-Step Review

**Date:** 2026-02-17
**Scope:** The full analysis path from raw recording to each primary output figure, covering assumptions, data transformations, intermediate storage, and known issues at every step.

---

## Table of Contents

1. [Entry Point and Raw Data](#1-entry-point-and-raw-data)
2. [Data Loading — `load_protocol2_data`](#2-data-loading)
3. [Bar Sweep Pipeline — Direction Selectivity](#3-bar-sweep-pipeline)
4. [Flash Pipeline — Receptive Field Mapping](#4-flash-pipeline)
5. [Bar Flash Pipeline — Spatial RF along PD/ND Axis](#5-bar-flash-pipeline)
6. [Newer Single-Experiment Pipeline — `analyze_single_experiment`](#6-newer-single-experiment-pipeline)
7. [Batch Pipeline — `batch_analyze_1DRF`](#7-batch-pipeline)
8. [Cross-Cutting Concerns](#8-cross-cutting-concerns)

---

## 1. Entry Point and Raw Data

### Entry point: `process_protocol2()`

**How it is called:** User `cd`s into an experiment folder and runs the function. There are no arguments — everything is inferred from the working directory.

```
process_protocol2()
  ├── loads currentExp.mat → metadata struct {Frame, Age, Strain, Side}
  ├── process_bars_p2(exp_folder, metadata, PROJECT_ROOT)   → resultant_angle
  ├── process_flash_p2(exp_folder, metadata, PROJECT_ROOT, resultant_angle)
  └── process_bar_flashes_p2(exp_folder, metadata, PROJECT_ROOT)
```

### Raw data structure

The G4 arena records a single continuous binary stream saved as `G4_TDMS_Log*.mat`. Inside:

| Field | Content | Units |
|-------|---------|-------|
| `Log.ADC.Volts(1,:)` | Frame position signal | Arbitrary integers — encodes which stimulus frame is being displayed |
| `Log.ADC.Volts(2,:)` | Membrane voltage | Volts (raw); multiplied by ×10 everywhere in the analysis |

**Sampling rate:** 10 kHz (1 sample = 0.1 ms). This is hard-coded throughout.

**Assumption:** The ×10 scaling converts the amplifier's output range to millivolts. This factor is never defined as a constant — it is a literal `*10` in every file that touches voltage.

### Stimulus order within one repetition (×3 reps total)

```
10s grey screen
4-pixel flashes (196 flashes, 14×14 grid)
6-pixel flashes (100 flashes, 10×10 grid)
3s grey screen
Bar sweeps — 28 dps  (16 directions)
3s grey screen
Bar sweeps — 56 dps  (16 directions)
3s grey screen
Bar sweeps — 168 dps (16 directions)
3s grey screen
Bar flashes — 80 ms  (8 orientations × 11 positions)
3s grey screen
Bar flashes — 14 ms  (8 orientations × 11 positions)
```

This sequence is the **fundamental structural assumption** of every parser. If the protocol changes, all parsers will break silently.

---

## 2. Data Loading

### `load_protocol2_data(exp_folder)`

**What it loads:**

| Source file | Variable loaded | Purpose |
|-------------|-----------------|---------|
| `Log Files/G4_TDMS_Log*.mat` | `Log` | Raw voltage + frame traces |
| `params/*.mat` (first file found) | `params` | Stimulus parameters: `.x`, `.y`, `.on_off`, crop sizes |
| `Functions/0001*.mat` | `pfnparam` | Position function metadata |

**Extracts from folder name:**
- `date_str` = characters `[end-15 : end-6]` of the folder name → `YYYY_MM_DD`
- `time_str` = characters `[end-4 : end]` → `HH_MM`

**Assumptions:**
- Folder name always ends with `YYYY_MM_DD_HH_MM` (16 characters)
- Path separator detection uses `contains(path, '/')` vs `'\'` — fragile on mixed-separator systems
- `params` always loads the **first** `.mat` file in the `params/` folder (OS-dependent sort order). Comment in code acknowledges this only gets 4px params.

**Side effect:** Changes working directory via `cd()` four times. Not restored on error. The newer `analyze_single_experiment` mitigates this with `onCleanup`.

**Stored/passed:** Returns `[date_str, time_str, Log, params, pfnparam]` — nothing saved to disk at this stage.

---

## 3. Bar Sweep Pipeline — Direction Selectivity

### Overview

```
process_bars_p2
  ├── load_protocol2_data  → Log, params
  ├── parse_bar_data       → data (32×4 cell)
  ├── plot_timeseries_polar_bars  → [max_v, min_v] + Figure 1
  ├── plot_polar_with_arrow       → resultant_angle + Figure 2
  ├── align_data_by_seq_angles    → data_ordered
  ├── find_PD_and_order_idx (×3 speeds) → d, ord, magnitude, angle_rad, fwhm, cv, thetahat, kappa
  ├── compute_bar_response_metrics (×3)  → sym_ratio, DSI, DSI_pdnd, vector_sum
  └── save → bar_results_*.mat
```

### Step 3a: `parse_bar_data(f_data, v_data)`

**Purpose:** Segment the continuous recording into individual bar sweep responses.

**Algorithm:**
1. Find all runs of `f_data == 0` (grey screen intervals)
2. Keep only runs ≥ 30,000 samples (≥ 3 seconds at 10 kHz)
3. These "3s gaps" are landmarks. The code indexes into them by ordinal position:
   - `idx_3(2)` to `idx_3(5)` = rep 1 bar sweeps
   - `idx_3(8)` to `idx_3(11)` = rep 2
   - `idx_3(14)` to `idx_3(17)` = rep 3

**Assumption:** Exactly 17+ zero-runs of ≥ 3s exist with this specific spacing. No validation.

4. Within each rep range, find transitions where `abs(diff(frames)) > 9` → these are stimulus onset/offset boundaries
5. Extract every-other segment (moving bar only, not inter-stimulus intervals), with ±9000 sample (0.9s) flanking windows
6. Trim all 3 reps to the minimum length, then compute `nanmean` across reps

**Output:** `data` = cell array, rows = 48 bar conditions (16 directions × 3 speeds), columns = [rep1, rep2, rep3, mean]

**Assumptions & risks:**
- The ±9000 sample flanking window is hard-coded. If stimuli are shorter than expected, flanks could overlap adjacent stimuli.
- `abs(diff) > 9` threshold for detecting transitions assumes frame numbers always jump by >9 between stimuli. A protocol with more frames could break this.

### Step 3b: `plot_timeseries_polar_bars` → **Figure 1: Radial timeseries + polar plot**

**What it does:**
1. Arranges 16 small timeseries plots in a circle using manual `axes('Position', ...)` calls
2. Reorders data rows via `plot_order = [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16]` to map from forward/backward pairs to sequential angular order
3. For each direction, plots individual reps in grey, mean in colour
4. Extracts `max_v` and `min_v` from the mean trace:

**Peak extraction (per direction, per speed):**
```
d_before_flash = mean_trace(1000:9000)        % 100–900 ms before stimulus
mean_before    = mean(d_before_flash)           % baseline
d_stim         = mean_trace(9000:end-7000)      % stimulus period, trimming 700ms at end
max_v(i, sp)   = abs(prctile(d_stim, 98) - mean_before)
min_v(i, sp)   = prctile(d_stim(half:end), 2)  % min in 2nd half only
```

**Assumptions:**
- Baseline is always in samples 1000–9000 (100–900 ms from trace start)
- 98th percentile is used instead of absolute max to reject noise spikes
- Min is taken only from the 2nd half of the stimulus period to avoid onset transients
- y-limits hard-coded to `[-80, -10]` — if a cell has unusual resting potential, traces may be clipped visually (does not affect data)

5. Central polar plot: `polarplot(angles, max_v_polar)` for each speed

**Stored:** `max_v` (16×3) and `min_v` (16×3) returned to caller
**Saved:** PDF to `figures/bar_stimuli/baselineb4/<on_off>/`

### Step 3c: `plot_polar_with_arrow` → **Figure 2: Polar DS plot with vector sum arrow**

For each of 3 speeds:
1. Subtracts `median_voltage` from max_v_polar (note: **not** done in Figure 1's polar plot — these use raw max_v)
2. Calls `vector_sum_polar` → computes resultant angle
3. Calls `add_arrow_to_polarplot` → overlays grey arrow

**The last speed's resultant_angle is returned** as the function output and passed to the flash pipeline.

**Assumption/risk:** The `resultant_angle` returned is always from the **3rd (168 dps) speed**, not the slow speed. But `find_PD_and_order_idx` in the next step uses the slow-speed data. The two could disagree about PD.

**Saved:** PNG + PDF to `figures/bar_stimuli/baselineb4/<on_off>/`

### Step 3d: `find_PD_and_order_idx(max_v_polar, median_voltage)` — called 3× (once per speed)

**Vector sum calculation:**
```matlab
angles    = linspace(0, 2π, 17)  →  [0, π/8, π/4, ... 2π]  (remove last)
responses = max_v_polar(1:16)
x_comp    = responses .* cos(angles)
y_comp    = responses .* sin(angles)
magnitude = sqrt(Σx² + Σy²) / Σ(responses)       % normalized 0–1
angle_rad = atan2(Σy, Σx)                          % preferred direction
```

**Assumption:** Responses are **not** baseline-subtracted in the vector sum here (the `- median_voltage` is commented out on line 47). This means the vector sum is computing with absolute response amplitudes, not relative depolarizations. This is inconsistent with Figure 2's arrow, which *does* subtract median.

**Tuning width metrics:**
- `compute_FWHM`: Finds all angles where response ≥ half-max, returns angular width between first and last. **Risk:** If the tuning curve has two peaks (bimodal), this will report the span between them, not the width of each peak.
- `compute_circular_var`: `CV = 1 - |vector_sum| / Σ(responses)`. Standard formula.
- `circ_vmpar`: Von Mises fit from the Circular Statistics Toolbox (external dependency).

**PD alignment:**
- Finds the discrete direction closest to the vector sum angle
- Computes `ord` — a 16-element reordering that maps PD to subplot position 8 (i.e., ~π/2)
- Creates `d` = 16×2 matrix `[aligned_angles, responses]`, sorted by angle

**Stored:** `d`, `ord`, `magnitude`, `angle_rad`, `fwhm`, `cv`, `thetahat`, `kappa` — all returned to caller

### Step 3e: `compute_bar_response_metrics(d)` — called 3× (once per speed)

**Inputs:** `d` = 16×2 [aligned_angles, responses], with PD at row 5 (π/2).

**Symmetry ratio:**
```
Compare 7 pairs of responses symmetric about PD:
  (row 4, row 6), (3,7), (2,8), (1,9), (16,10), (15,11), (14,12)
sym_val   = Σ|diff| / Σ(all responses)
sym_ratio = 1 - sym_val          % 1 = perfectly symmetric
```

**DSI (vector sum):**
```
vector_sum = Σ(response × exp(i × angle))
DSI        = |vector_sum| / Σ(responses)
```

**DSI (PD/ND):**
```
DSI_pdnd = (d(5,2) - d(13,2)) / (d(5,2) + d(13,2))
```
where row 5 = PD and row 13 = ND (opposite direction, 180° away).

**Assumption:** Row 5 is always PD and row 13 is always ND. This is guaranteed by the alignment step in `find_PD_and_order_idx`, but only if the alignment worked correctly.

### Step 3f: Saving bar results

**File:** `bar_results/<strain>_<on_off>_<date>_<time>.mat`

**Contents:**

| Variable | Shape | Description |
|----------|-------|-------------|
| `bar_results` | struct | All metrics: Date, Time, Strain, Type, median_voltage, and for each speed: max_v_polar, magnitude, angle_rad, fwhm, cv, thetahat, kappa, sym_ratio, vector_sum, DSI_vector, DSI_pdnd |
| `data` | 48×4 cell | Raw parsed bar data (all reps + mean, all speeds) |
| `data_aligned` | 48×4 cell | PD-aligned version |
| `ord` | 16×1 | Reordering indices |
| `d_slow`, `d_fast`, `d_vfast` | 16×2 each | PD-aligned [angles, responses] per speed |

---

## 4. Flash Pipeline — Receptive Field Mapping

### Overview

```
process_flash_p2
  ├── load_protocol2_data → Log, params
  ├── Find flash start indices (different logic for ON vs OFF)
  ├── For px_size = [4, 6]:
  │     ├── parse_flash_data → data_comb, cmap_id, var metrics, max_data, min_data
  │     ├── rescale(data_comb, 0, 1) → data_comb2
  │     ├── plot_rf_estimate_timeseries_line → Figure 3 (spatial timeseries grid)
  │     ├── plot_heatmap_flash_responses     → Figure 4 (RF heatmap)
  │     ├── gaussian_RF_estimate → optEx, R², optInh, R²_inh + Figures 5, 6
  │     └── store in rf_results struct
  └── save → rf_results_*.mat
```

### Step 4a: Finding flash start indices

**ON vs OFF flashes use completely different detection logic:**

**OFF flashes:**
```matlab
idx = find(diff_f_data == 1 & f_data(2:end) == 1)  % frame jumps to 1
idx = idx([1,2,5,6,9,10])  % pick specific indices
```

**ON flashes:**
```matlab
idx_4 = find(diff_f_data == 197 & f_data == 197)  % 4px: frame jumps to 197
idx_6 = find(diff_f_data == 101 & f_data == 101)  % 6px: frame jumps to 101
idx = sort([idx_4, idx_6])
```

**Assumption:** Hard-coded frame values (1, 197, 101) depend on the exact number of flash patterns. If the protocol adds or removes patterns, these will silently select wrong indices.

### Step 4b: `parse_flash_data` — per-flash response extraction

**For each of the 196 (or 100) flashes, for each of 3 reps:**

1. Find all positive frame transitions within the rep's range → `start_flash_idxs`
2. Extract a 7000-sample (0.7s) window: 1000 samples before flash onset + 6000 after
3. Compute mean across 3 reps
4. Classify the response:

**Response classification:**
```
max_val = prctile(mean_response(500:end), 98)       % peak in excitatory direction
min_val = prctile(mean_response(2500:end), 2)        % trough in inhibitory direction
diff_resp = max_val - min_val

if |max| ≥ |min| AND diff_resp > 3:
    val = max_val,  type = excitatory (1)
elif |max| < |min| AND diff_resp > 2.8:
    val = min_val,  type = inhibitory (2)
else:
    val = mean(last 25% of trace),  type = neutral (3)
```

**Assumptions:**
- Excitatory threshold = 3, inhibitory threshold = 2.8 (in mV after ×10 scaling). These are **hard-coded magic numbers** with no documentation of how they were chosen.
- The asymmetric thresholds (3 vs 2.8) suggest inhibitory responses are typically smaller; this may not hold for all cell types.
- For slow flashes, max is searched from sample 500 onwards (50ms after window start = flash onset). Min is searched from sample 2500 (250ms). This assumes the inhibitory response is always delayed relative to excitation.
- The TODO on line 218 indicates the baseline period (first 1000 samples) was added later and the thresholds may not have been updated.

5. **Grid position mapping** (frame number → row, col):
```matlab
% OFF flashes (frames 0–195):
rows = 14 - mod(frame_num, 14)
cols = floor(frame_num / 14) + 1

% ON flashes (frames 196–391):
rows = 14 - mod(frame_num - 196, 14)
cols = floor((frame_num - 196) / 14) + 1
```

**Stored per flash position:** `data_comb(row,col)` = response value, `cmap_id` = classification, variance metrics, max/min data.

### Step 4c: Rescaling

```matlab
data_comb2 = rescale(data_comb, 0, 1)
```

The entire NxN response matrix is linearly rescaled to [0, 1]. This is used for the heatmap and Gaussian fitting. **The original data_comb (in mV) is preserved separately in the results struct.**

### Step 4d: `plot_rf_estimate_timeseries_line` → **Figure 3: Spatial timeseries grid**

Creates a 14×14 (or 10×10) grid of tiny subplots. Each subplot shows the averaged voltage trace for that flash position. Red colouring indicates excitatory response strength. An arrow overlay from the bar analysis indicates the cell's preferred direction.

**Saved:** PDF + PNG

### Step 4e: `plot_heatmap_flash_responses` → **Figure 4: RF heatmap**

Simple `imagesc(data_comb2)` with a red-blue diverging colormap (`redblue`). Color limits centered on median ± 0.5.

**Saved:** PDF

### Step 4f: `gaussian_RF_estimate` → **Figures 5 & 6: Gaussian RF fits**

**Excitatory lobe fitting:**
1. Create coordinate grid matching the NxN response matrix
2. Log-transform the rescaled data: `z = sign(z) × log(1 + |z|)`
3. Define a rotated 2D Gaussian with 7 parameters: `[A, x₀, y₀, σx, σy, θ, B]`
4. Fit using `lsqcurvefit` with bounds

**Inhibitory lobe fitting:**
1. Takes the original (non-rescaled) `min_data` matrix, inverts it (`× -1`)
2. Fits the same Gaussian model (but does NOT log-transform — uses raw inverted values)
3. **Note:** The excitatory fit uses `response` (rescaled 0–1, then log-transformed), while the inhibitory fit uses `min_data × -1` (raw mV). This asymmetry is potentially intentional (different dynamic ranges) but undocumented.

**Figure 5:** Side-by-side comparison: original data | excitatory fit | inhibitory fit
**Figure 6:** RF heatmap with 1.5σ contour ellipses overlaid (red = excitatory, black = inhibitory)

**Note:** Figure 6 calls `close` at line 171, which closes the figure immediately after creating it. The figure handle `f2` is returned but the figure no longer exists visually.

**Assumption:** The Gaussian initial guess uses `mean(xData)` for the center, which is the grid center. If the RF is near an edge, the optimiser may converge to a local minimum. No multi-start fitting or robustness check is performed.

### Step 4g: Saving flash results

**File:** `rf_results_<date>_<time>_<strain>_<on_off>.mat`

**Contents (per speed — currently only "slow"):**

| Field | Description |
|-------|-------------|
| `.data_comb` | NxN response values (mV) |
| `.max_data`, `.min_data` | NxN peak/trough values |
| `.cmap_id` | NxN response classification |
| `.var_within_reps`, `.var_across_reps` | NxN CV matrices |
| `.var_filtered_v` | Variance of moving-mean-filtered voltage (overall quality metric) |
| `.R_squared`, `.R_squaredi` | Gaussian fit quality |
| `.sigma_x_exc`, `.sigma_y_exc` | Excitatory RF extent (pixels) |
| `.sigma_x_inh`, `.sigma_y_inh` | Inhibitory RF extent (pixels) |
| `.optExc`, `.optInh` | Full 7-parameter Gaussian fit vectors |

---

## 5. Bar Flash Pipeline — Spatial RF along PD/ND Axis

### Overview

```
process_bar_flashes_p2
  ├── load_protocol2_data → Log, params
  ├── parse_bar_flash_data → data_slow, data_fast, mean_slow, mean_fast
  ├── plot_bar_flash_data (80ms) → Figure 7
  ├── plot_bar_flash_data (14ms) → Figure 8
  └── save → bar_flash_results_*.mat
```

### Step 5a: `parse_bar_flash_data(f_data, v_data)`

**Structure:** 8 bar orientations × 11 spatial positions × 2 speeds × 3 reps

**Algorithm:**
1. Find ≥3s zero-runs (same as bar sweep parser)
2. For slow bar flashes: rep 1 = between `idx_3(5)` and `idx_3(6)`, rep 2 = `idx_3(11:12)`, rep 3 = `idx_3(17:18)`. For fast: offset by 1.
3. Within each rep, find flash onsets (`diff(f_data > 0)`) and pair starts/ends
4. Use `max(f_data(start:end))` as the flash frame number → determines position in the `data{position, orientation}` cell
5. Extract voltage trace ± `gap_between_flashes` around each flash (5000 samples for slow = 0.5s, 2500 for fast)
6. Compute mean across reps with NaN-padding for unequal lengths

**Output:**
- `data_slow`, `data_fast`: 11×8×3 cell arrays (positions × orientations × reps)
- `mean_slow`, `mean_fast`: 11×8 cell arrays (mean traces)

**Assumption:** `data_rep{flash_frame_num} = data_flash` — this directly indexes the cell array using the frame number. If frame numbers don't map sequentially from 1 to 88 for each speed, data will be placed in wrong cells or cause index-out-of-bounds.

### Step 5b: `plot_bar_flash_data` → **Figures 7 & 8: Bar flash grids**

Creates an 8×11 tiled layout:
1. Computes background colour from the 98th percentile of the mean trace (50%–75% of trace duration) per position, normalized to max across all positions
2. Plots individual reps in grey, mean in black
3. Red background intensity indicates response magnitude

**Saved:** PDF

### Step 5c: Saving bar flash results

**File:** `bar_flash_results_<date>_<time>_<strain>_<on_off>.mat`

**Contents:** `data_slow`, `data_fast`, `mean_slow`, `mean_fast` (raw cell arrays)

---

## 6. Newer Single-Experiment Pipeline — `analyze_single_experiment`

This is a **cleaner rewrite** of parts of the above pipeline with several improvements:

**Key differences from `process_protocol2`:**
1. Uses `opts` struct for all configurable parameters (with defaults)
2. Loads `bar_lut.mat` — a lookup table mapping pattern IDs to directions/orientations
3. Uses `onCleanup` to restore working directory
4. Validates LUT directions against the actual `currentExp.mat` pattern/function order
5. Uses `compute_bar_sweep_responses` (baseline-subtracted depolarization) instead of inline percentile extraction
6. Uses `find_pd_from_lut` (LUT-aware PD finding with explicit orthogonal mapping) instead of `find_PD_and_order_idx`
7. Generates bar flash figures along both PD-ND and orthogonal axes

**Figure outputs:**
- Slow bar sweep polar timeseries with vector sum arrow
- 8×11 bar flash heatmap (all orientations)
- 1×11 bar flash subplots along PD-ND axis
- 1×11 bar flash subplots along orthogonal axis

**All saved as 300 dpi PDF** to `<exp_folder>/analysis_output/`

---

## 7. Batch Pipeline — `batch_analyze_1DRF`

Processes all experiment folders in a data root directory:

1. Discovers experiment folders (date-formatted subfolders)
2. For each experiment: runs `analyze_single_experiment` logic, extracts baseline-subtracted bar flash traces
3. Classifies cells as ON/OFF (based on `metadata.Frame > on_threshold` where threshold = 129) and control/ttl (based on strain name containing 'ttl')
4. Groups cells into 4 categories: on_control, on_ttl, off_control, off_ttl
5. Generates population-averaged plots (3 per ON/OFF condition = 6 total):
   - PD-aligned polar tuning curves (control = black, ttl = red, ± spread)
   - 1×11 PD-ND bar flash traces
   - 1×11 orthogonal bar flash traces
6. Spread = MAD (default) or SEM, configurable via `opts.stat_method`

**Output:** `results` struct array with per-cell data + population PDFs

---

## 8. Cross-Cutting Concerns

### 8.1 Baseline handling inconsistencies

| Pipeline step | Baseline definition | Method |
|---------------|-------------------|--------|
| `plot_timeseries_polar_bars` | samples 1000–9000 | `mean()` of pre-stimulus window |
| `find_PD_and_order_idx` | None applied | Vector sum uses raw `max_v_polar` (no baseline subtraction) |
| `plot_polar_with_arrow` | `median_voltage` | Subtracts recording-wide median before vector sum |
| `compute_bar_sweep_responses` | samples 1000–9000 | `mean()` of per-direction baseline |
| `parse_flash_data` | `median(v_data)` globally | Subtracts recording-wide median, then classifies with absolute thresholds |

The **recording-wide median** vs **per-stimulus baseline** difference means metrics from different parts of the pipeline are not directly comparable. If membrane potential drifts during a long recording, the global median will be different from local baselines.

### 8.2 `max_v` definition varies between old and new pipelines

| Pipeline | max_v calculation |
|----------|------------------|
| `plot_timeseries_polar_bars` (old) | `abs(prctile(stim_period, 98) - mean(baseline))` |
| `compute_bar_sweep_responses` (new) | `abs(prctile(stim_period, 98) - mean(baseline))` — same formula but with configurable ranges |

These are consistent in formula but differ in hard-coded vs configurable sample ranges.

### 8.3 The `plot_order` assumption

```matlab
plot_order = [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16]
```

This mapping is defined identically in 4 separate files (`plot_timeseries_polar_bars.m`, `align_data_by_seq_angles.m`, `analyze_single_experiment.m`, and comments elsewhere). It reflects how the G4 arena presents bars in forward/backward pairs. If the stimulus presentation order ever changes, all 4 files must be updated.

### 8.4 Response classification thresholds

The flash classification thresholds (`diff_resp > 3` for excitatory, `> 2.8` for inhibitory) are hard-coded in `parse_flash_data.m` with no documented justification. These operate on median-subtracted, ×10-scaled voltage, so they correspond to ~0.3 mV and ~0.28 mV real membrane potential changes. Changes to amplifier gain or recording conditions could make these thresholds inappropriate.

### 8.5 Figure close/save inconsistencies

- `gaussian_RF_estimate` calls `close` on figure `f2` immediately after creating it (line 171). The figure handle is returned but the figure window is gone.
- `plot_polar_with_arrow` saves PNG/PDF to the **current directory** (not `fig_folder`), despite constructing `fig_folder` and creating it if missing. The `fig_folder` is never used in the `exportgraphics` calls.
- Some figures are saved as vector PDF, others as raster PNG, and some as both. No consistent policy.

### 8.6 Data flow summary

```
Raw TDMS recording
       │
       ▼
 load_protocol2_data ──→ Log.ADC.Volts, params, pfnparam
       │
       ├──→ parse_bar_data ──→ 48×4 cell (voltage traces per direction × speed)
       │         │
       │         ├──→ plot_timeseries_polar_bars  ──→ max_v (16×3), min_v (16×3)
       │         │         │                            │
       │         │         └──→ Figure 1: Radial timeseries + polar
       │         │
       │         ├──→ plot_polar_with_arrow ──→ resultant_angle
       │         │         │
       │         │         └──→ Figure 2: Polar DS + arrows
       │         │
       │         ├──→ find_PD_and_order_idx ──→ d (16×2), metrics
       │         │
       │         └──→ compute_bar_response_metrics ──→ DSI, sym_ratio
       │
       │         ═══ SAVED: bar_results_*.mat ═══
       │
       ├──→ parse_flash_data (×2: 4px, 6px) ──→ NxN response grids
       │         │
       │         ├──→ plot_rf_estimate_timeseries_line
       │         │         └──→ Figure 3: Spatial timeseries grid
       │         │
       │         ├──→ plot_heatmap_flash_responses
       │         │         └──→ Figure 4: RF heatmap
       │         │
       │         └──→ gaussian_RF_estimate
       │                   ├──→ Figure 5: Gaussian fits comparison
       │                   └──→ Figure 6: RF + 1.5σ contours
       │
       │         ═══ SAVED: rf_results_*.mat ═══
       │
       └──→ parse_bar_flash_data ──→ 11×8×3 cell arrays
                 │
                 ├──→ plot_bar_flash_data (80ms)
                 │         └──→ Figure 7: 8×11 bar flash grid (slow)
                 │
                 └──→ plot_bar_flash_data (14ms)
                           └──→ Figure 8: 8×11 bar flash grid (fast)

               ═══ SAVED: bar_flash_results_*.mat ═══
```
