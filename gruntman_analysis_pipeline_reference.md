# Gruntman et al. (Reiser Lab) — Analysis Pipeline Reference

> Extracted from: Gruntman, Romani & Reiser, *Nat Neurosci* 2018; *eLife* 2019; Gruntman, Reimers, Romani & Reiser, *Curr Biol* 2021.
> For use as context in a Claude Code session for replicating or extending these analyses.

---

## Source Code & Data Links

| Resource | URL |
|----------|-----|
| T5 conductance model (MATLAB + all 17 cell .mat files) | https://github.com/reiserlab/T5ConductanceModel |
| eLife 2019 figure data + plotting code (figshare) | https://doi.org/10.25378/janelia.c.4771805.v1 |
| Nat Neurosci 2018 data | https://doi.org/10.25378/janelia.5576101 |
| Nat Neurosci 2018 analysis code (figshare) | https://figshare.com/projects/Gruntman_et_al_2017_Data/26347 |

---

## 1. Recording & Acquisition Parameters

- **In-vivo whole-cell current clamp**, Axon MultiClamp 700B
- **Sampling**: 20 kHz, hardware low-pass at 10 kHz
- **Holding potential**: ~−65 mV (uncorrected for liquid junction potential), maintained with 0–3 pA constant hyperpolarizing current
- **Cell types**: T4 (ON, split-GAL4 SS02344) and T5 (OFF, split-GAL4 SS25175)
- **Display**: Modular LED panels, 216° × 72°, pixel ≈ 2.25°, green LEDs (565 nm), background ~31 cd/m², bright ~72 cd/m²
- **Cell exclusion**: Recording terminated if visual or current-step responses diminished noticeably (no numerical threshold)

### Undocumented software filter (from figshare code)

`fdatooLowPass.m` defines a 60 Hz passband / 500 Hz stopband equiripple FIR filter at 10 kHz Fs. **Not mentioned in any paper.** If applied to analysis data (not just display), it would limit temporal resolution to ~16 ms. Context of application unknown — main analysis pipeline was not included in figshare uploads.

---

## 2. Baseline Subtraction

- Each trace has **900 ms of pre-stimulus baseline** (18,000 samples at 20 kHz; `newZero = floor(45000/2.5)`)
- Data in figshare `.mat` files is **already baseline-subtracted** — stored as mV deviation from rest (y-range: −7.5 to +25 mV), not absolute Vm around −65 mV
- In the conductance model: `vm_output = vm - Vleak` where `Vleak = -65 mV`
- **Exact baseline window** (e.g., mean of how many ms before onset) is not specified in papers or available code — computed in unpublished preprocessing scripts upstream of the figure-plotting code

---

## 3. Trial Exclusion

- **Papers 1 & 2**: No explicit numerical criteria stated
- **Paper 3 (Curr Biol 2021)** specifies:
  - Exclude trial if pre-stimulus baseline mean differs from overall pre-stimulus mean for that stimulus group by > **10 mV**
  - Exclude trial if pre-stimulus mean and overall trial mean differ by > **15 mV** (or > **25 mV** for slow moving bars)
- Typically **n = 3 trials** per condition per cell
- Data is **trial-averaged per cell before** population averaging (figshare `.mat` stores one trace per cell per condition)

---

## 4. RF Center Localization

Hierarchical, stepwise procedure:

1. Present coarse grid of OFF (T5) or ON (T4) flashing squares (~11° × 11°) spanning ~79° × 68°
2. Identify grid position with largest depolarization
3. Present smaller stimuli at higher resolution within responsive region
4. Repeat until peak localized to **single pixel (~2° × 2°)** → this becomes **position 0**

---

## 5. PD-ND Axis Identification

- Present narrow bar (1 × 9 LEDs, ~2.25° × 20.25°) moving in **8 directions** at 28°/s through RF center
- Direction with largest depolarization = **Preferred Direction (PD)**
- Opposite = **Null Direction (ND)**
- All subsequent mapping along the identified PD-ND axis

---

## 6. Single Position Flash Responses (SPFRs) — First-Order RF

### Stimulus
- Bars of width 1, 2, 4 pixels (2.25°, 4.5°, 9°) × 9 pixels long
- Flash durations: **40 ms** (short) and **160 ms** (long) — Papers 2 & 3; Paper 1 also uses 20, 80, 320 ms
- OFF bars for T5 (pixels OFF from background); ON bars for T4 (pixels ON from background)
- Non-preferred contrast (NC) bars in Paper 3: OFF for T4, ON for T5
- Positions randomly interleaved along PD-ND axis

### Alignment across cells
- All cells aligned to **position of peak depolarization = position 0**
- Negative positions = leading side (PD enters from here)
- Positive positions = trailing side (ND enters / hyperpolarization detected)
- 1 position unit = 1 pixel ≈ 2.25° (exact angle depends on cardinal vs diagonal PD — see eLife 2019 Fig 2—supplement 1)

### Position correction for bar width (from code)
```matlab
corrPInd = pInd + floor(width / 1.5);
% width 1 → +0; width 2 → +1; width 4 → +2
```

### Population averaging
```matlab
plotMean = nanmean(relDat, 2);    % mean across cells
plotSEM  = nanstd(relRedDat, 0, 2) / sqrt(size(relRedDat, 2));  % SEM across cells
```
- `nanmean`/`nanstd` handles cells missing certain positions
- SEM reflects **between-cell variability** (within-cell trials already averaged)

### Edge trimming
- 9 padded positions trimmed from each edge: `chopEdges = 9`
- Leaves ~13 plotted positions spanning the RF

---

## 7. RF Metrics Extracted from SPFRs

### Peak responses (per cell, per position)
- **Peak depolarization**: `max(trace)` during stimulus + post-stimulus window
- **Peak hyperpolarization**: `min(trace)` during same window
- Population summary: **median** with **upper and lower quartiles** (eLife 2019 Fig 2D)

### Temporal metrics (eLife 2019, from 160 ms width-2 bar flashes)
- **Rise start time**: time to reach **10% of max** depolarization
- **Rise time**: time from **10% to 50%** of max
- **Decay time**: time from **80% to 20%** of max
- Reported as **difference from central position** (position 0)
- Stats: two-tailed t-test, p < 0.05, pooled across leading/trailing sides

### Onset time (Nat Neurosci 2018)
- Constant across RF (no evidence for HR-type temporal delay)
- Significant position dependence in **decay time** (faster on trailing side due to inhibition)
- Stats: slope of linear regression of time vs position; one-sided unpaired t-test, p < 0.01

## 8. Moving Bar Responses & DSI

### Stimulus construction
- Moving bar = ordered sequence of single-position bar flashes
- Bar widths: 1, 2, 4 pixels; lengths: 9 pixels
- Speeds: 14°/s (160 ms/pixel), **28°/s (80 ms/pixel)** [primary], 56°/s (40 ms/pixel), 112°/s (20 ms/pixel)
- 8 directions through RF center; all cells aligned to their individual PD

### DSI formula
```
DSI = (PDmax − NDmax) / PDmax
```
- **PDmax**: max baseline-subtracted depolarization for PD motion
- **NDmax**: max baseline-subtracted depolarization for ND motion
- Both are peak Vm above baseline during the full stimulus presentation window

### Moving bar sample sizes (eLife 2019)
- n = 9 cells (width 1), 15 cells (width 2), 14 cells (width 4) — not all 17 cells received all widths
- Stats (Nat Neurosci 2018): DSI > 0, one-sided unpaired t-test, p < 0.05

### Linear prediction comparison
```
Linear_PD = Σ(SPFRs aligned in PD temporal order)
Linear_ND = Σ(SPFRs aligned in ND temporal order)
```
- The only difference between summed PD and ND is the temporal order of the same SPFRs
- Linear model consistently **overestimates** measured responses (especially ND and wide/fast bars)

---

## 9. Two-Step Apparent Motion (Second-Order RF)

### Stimulus
- Two bars flashed sequentially at adjacent or overlapping positions
- Same widths (2, 4) and durations (40, 160 ms) as SPFR mapping
- PD sequence: trailing → leading position; ND sequence: reversed

### Analysis
Compare **measured** two-bar response vs **summed** (superposition of two individual SPFRs, time-aligned):

| Metric | Computation | Indicates |
|--------|-------------|-----------|
| PD enhancement | Measured_PD − Summed_PD > 0 | Amplification in PD |
| ND suppression | Measured_ND − Summed_ND < 0 | Suppression in ND |

Two comparison metrics used (both give same conclusion):
1. **Response maximum** (peak during 2nd bar window)
2. **Response mean** (mean during 2nd bar window) — shown in Fig 3—supplement 1

### Key result
- **No PD enhancement** found under any condition (widths 2 or 4, durations 40 or 160 ms, adjacent or overlapping)
- **ND suppression** correlates with DS magnitude (R² = 0.40, slope CI: [−0.64, −0.31])
- Minimum **n ≥ 3 cells** required per position pair for inclusion

### Statistics
- Unpaired t-test corrected for multiple comparisons via **FDR with q = 0.075**
- Linear regression with **95% CI on slope**

---

## 10. Grating Stimuli

### Static grating flashes (eLife 2019)
- Square wave (dark + background); multiple phases; durations 40 and 160 ms
- **Only cardinal-PD cells** used (n = 5 of 17) — grating must align with display grid

### Drifting gratings
- Two speeds: TF = 3.125 Hz (40 ms steps) and 0.78 Hz (160 ms steps)
- Two starting phases; PD and ND directions
- Response characterized by **fit amplitude** and **fit phase** (fitting method unspecified — likely sinusoidal least-squares)

### Cells with grating data (from model code, line 36)
- Cells [2, 5, 6, 7, 9, 10, 13, 14, 15, 16, 17] — 11 of 17 total

---

## 11. Conductance Model (EI Model) — Full Specification

### Membrane equation
```matlab
Vm = (Vl + Ve*ge + Vi*gi) / (1 + ge + gi)
% Output: Vm - Vl (deviation from rest)
```

### Fixed parameters
| Parameter | Value |
|-----------|-------|
| V_leak | −65 mV |
| V_exc | 0 mV |
| V_inh | −74 mV |
| Pre-stim offset | 25 ms |

### Temporal dynamics (ODE system)
```matlab
dhe/dt = (-he + ae) / tre     % excitatory hidden state
dge/dt = (-ge + he) / tde     % excitatory conductance (output)
dhi/dt = (-hi + ai) / tri     % inhibitory hidden state
dgi/dt = (-gi + hi) / tdi     % inhibitory conductance (output)
```
Cascaded first-order → alpha-function profile; 2 time constants per conductance.

### Spatial filters (Gaussians)
```matlab
xae = ae * exp(-(xx - mue)^2 / (2*sige^2))
xai = ai * exp(-(xx - mui)^2 / (2*sigi^2))
```

### 11 free parameters (with optimization bounds)

| # | Param | Description | Lower | Upper |
|---|-------|-------------|-------|-------|
| 1–5 | tau, tde, tdi, tre, tri (×10) | Time constants | 1 ms | 400 ms |
| 6–7 | ae, ai | Conductance amplitudes | 0 | 10 |
| 8 | mue | Excitatory center | −5 | 5 (pixels) |
| 9 | mui | Inhibitory center | −5 | **10** (pixels) |
| 10–11 | sige, sigi | Spatial sigmas | 0 | 10 |

Note: `mui` bound is asymmetric [−5, 10] — permits inhibitory field to be offset to trailing side.

### Optimization
- **Optimizer**: `fmincon` (MATLAB Optimization Toolbox)
- **Objective**: mean squared error between model and measured traces
- **Fit data**: width-2 bar flash responses ONLY (code line 53: `if(wid==2)`)
- **Random starts**: 10 per cell per seed (`attempts_num = 10`); paper reports running 1000 seeds on a cluster → top 10 solutions selected per cell
- **Early stop**: `resnorm < 3.5`
- **ODE solver**: `ode45`
- All other stimuli (width 1, width 4, moving bars, gratings) are **pure predictions**

### Prediction quality metric
- **MAD** (Mean Absolute Deviation) of measured responses to repeated stimuli
- ± upper quartile MAD plotted as bounds around identity line in predicted-vs-measured scatter plots

---

## 12. Key Discrepancies & Gaps

| Item | Details |
|------|---------|
| **Rearing conditions** | Papers 1 & 2: constant light. Paper 3: 16:8 L:D |
| **T5 driver line** | SS25175 named only in Paper 3; likely same line used in Paper 2 |
| **Trial exclusion** | Numerical thresholds (10/15/25 mV) stated only in Paper 3 |
| **Low-pass filter** | 60 Hz FIR in figshare code; **undocumented** in all papers |
| **Baseline window** | 900 ms pre-stimulus available; exact subtraction window unspecified |
| **Optimization** | Paper 1: all widths/speeds. Paper 2: width-2 only. Paper 3: composite ON-OFF |
| **Model complexity** | Paper 1: multi-compartment EM reconstruction. Paper 2: single-compartment. Paper 3: 4-conductance |
| **Paper 3 model code** | Not publicly released |
| **Grating fit method** | Amplitude/phase extraction procedure not detailed |
| **Rise time asymmetry** | T5 trailing-side rises faster (Paper 2); T4 has no onset difference (Paper 1) |
