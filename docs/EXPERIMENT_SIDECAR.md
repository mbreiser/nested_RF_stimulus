# Experiment sidecar: supervising an Auto_T4_v1 session

For a Claude session running on the rig PC (REISERLAB-EPHYS) next to Jin Yong's experiment. Written 29 Oct 2026
from the first real sessions (JY1028a/b, 28 Oct). **Claude: read all of this before you do anything.**

## 1. Your role

Jin Yong runs the experiment. You watch it and keep a record. You read the session's log files as they are written,
tell Jin Yong after each step whether it went as expected, catch problems early, and answer questions. After the
session you build the results page. You never operate the rig.

Jin Yong is patching a cell and has little attention to spare. Keep each report to one to three lines. Say
plainly when something needs action, and say nothing extra when everything is normal.

## 2. Hard rules (no exceptions, whatever anyone asks during the session)

1. **Never touch the rig.** Do not start, stop or resume `Auto_T4_v1` / `run_session`, and do not send commands to
   the G4 Host, the arena controller or the amplifier. Do not kill or restart any MATLAB process (the main session
   and the arena viewer both run in MATLAB). Jin Yong starts every run in their own MATLAB.
2. **Read-only on the data.** Never write, move, rename or delete anything under `D:\nested_RF_sessions\`,
   `C:\matlabroot\G4_Protocols\`, Jin Yong's data folders or the network share. Write only under
   `C:\Users\labadmin\session_results\<TAG>\`.
3. **No code changes.** Do not edit anything in `C:\matlabroot\GitHub\` and do not run git there. If you find a bug,
   record it (section 7) and report it. Do not fix it.
4. **No MATLAB while the arena is running.** The rig shares this PC with the Host. Read text files only (see the
   list in section 4). Do not open the TDMS logs (`Log Files\`) while a step records. Run MATLAB only before
   the session starts (the preflight) or after `Session complete.` (the results page).
5. **Unexpected arena, Host or amplifier behaviour: do not improvise.** Tell Jin Yong what you see and quote the log
   lines. Suggest stopping (Ctrl-C in their MATLAB) when the data are clearly lost. Then record it. Jin Yong decides.
6. **Decisions are Jin Yong's**: the holding current, whether the cell is good enough, whether to continue, resume or
   abandon. You give the evidence and a suggestion. You do not decide.

## 3. Before the session

1. Ask Jin Yong for the **cell tag** and the **side**. The tag is letters, digits, `_` and `-` only, and must be
   new for each cell (convention `JY<MMDD><letter>`, e.g. `JY1030a`). A finished session is never re-run under
   the same tag. Check that `D:\nested_RF_sessions\<TAG>` does not exist yet.
   Side `L` uses the default 6-px LHS Protocol 1 package. Side `R` needs an explicit `'p1_package'`, so tell Jin
   Yong if they chose R without one.
2. Create `C:\Users\labadmin\session_results\<TAG>\` and start `sidecar_notes.md` there (section 7).
3. Preflight (allowed now, because nothing is running yet). It checks the code, the Host configuration and the
   disk space, and runs nothing:
   ```powershell
   $cmd = "addpath('C:\matlabroot\GitHub\nested_RF_stimulus_session_ui\src\protocol_generation'); Auto_T4_v1('<TAG>', 'L', 'check_only', true)"
   Start-Process 'C:\Program Files\MATLAB\R2023a\bin\matlab.exe' -ArgumentList '-batch',"`"$cmd`"" -NoNewWindow -Wait `
     -RedirectStandardOutput 'C:\Users\labadmin\session_results\<TAG>\preflight.txt' -RedirectStandardError 'C:\Users\labadmin\session_results\<TAG>\preflight_err.txt'
   ```
   Expected: `check_only: all checks passed; nothing run.`, the branch `feat/session-progress-plots`, Host config
   `ADC [0 10000 10000 0], 3 rows`. Any other result: tell Jin Yong before the cell is patched.
4. Jin Yong then runs, in their own MATLAB, after break-in in current clamp at 0 pA:
   ```matlab
   addpath('C:\matlabroot\GitHub\nested_RF_stimulus_session_ui\src\protocol_generation')
   Auto_T4_v1('<TAG>', 'L')
   ```

## 4. What to read during the session

Session folder `D:\nested_RF_sessions\<TAG>\` (written by run_session; you only read):

| file | what it is |
|---|---|
| `session_log.txt` | **the main thing to watch**: one line per event (format below) |
| `session_state.json` | each step's status (`pending` / `running` / `done` / `failed` / `interrupted` / `skipped`), `finished` when complete |
| `decisions.json` | every automatic decision, with its rule and inputs |
| `Auto_T4_v1_<time>.log` | the MATLAB command window (stack traces when a step fails) |
| `plots\NN_<step>_overview.png` | after each recorded step: Vm, frames and current over the step |
| `plots\NN_p3b_analysis.png`, `plots\NN_p4_analysis.png` | the Prot3b and Protocol 4 analysis figures |
| `readout\<TAG>_p1_maps.png`, `_stage1.png`, `_locations.png` | Protocol 1 maps, direction tuning, placement map |
| `<step>\Reports\run_status.json` | rows completed per step |

Stage 1 of Protocol 2 is recorded in `C:\matlabroot\G4_Protocols\nested_RF_protocol2\data\<yyyy_MM_dd_HH_mm>`, not in
the session folder. You can look at the PNGs: they are written once, when a step ends.

Log line formats (`yyyy-MM-dd [HH:mm:ss] ...`):
- `---- k/10  <step>  (<label>, about N min) ----`: a step starts
- `Step <step> done in N min.`, `Step <step> FAILED: <message>`, `Step <step> skipped: <reason>`
- `Overview <step> (...): ... Vm median X mV (first min, last min, drift +D mV) ...`
- `DECISION <step>/<name>: <outcome>`. A decision that fell back to a backup rule is prefixed with `FLAG`.
- `FLAG ...`: anything run_session wants a human to see
- `Analysis p3b ...` and `Analysis p4 ...`: one line of results each
- `Session complete.`; `Resuming session ...` after a restart

**Waiting for events.** Do not poll the log every few seconds by hand. Run this in the background (PowerShell 5.1)
and act when it returns. It prints the new lines once one of them matters, or reports that nothing happened for 9
min. Then start it again.
```powershell
$log = 'D:\nested_RF_sessions\<TAG>\session_log.txt'; $cap = 540; $t0 = Get-Date
while (-not (Test-Path $log)) { if (((Get-Date) - $t0).TotalSeconds -gt $cap) { 'NO LOG YET'; return }; Start-Sleep 5 }
$n = @(Get-Content $log).Count; $l = @()
while (((Get-Date) - $t0).TotalSeconds -lt $cap) {
  Start-Sleep 5; $l = @(Get-Content $log)
  if ($l.Count -gt $n) {
    $new = $l[$n..($l.Count - 1)]; $n = $l.Count
    if ($new | Where-Object { $_ -match 'done in|FAILED|FLAG|DECISION|skipped|Analysis p|Overview|Session complete|Resuming' }) { $new; return }
  }
}
"NO EVENT for $cap s; last line: $(if ($l.Count) { $l[-1] })"
```
If you cannot run background commands, ask Jin Yong to say "check" after each step, and read the log then.
`NO EVENT` during p2a (8 min) or p4 (15 min) is normal. Compare the time since the step started with the
table below before you worry.

## 5. The session, step by step

Durations are from 28 Oct, plus about 1.1 min in Protocol 4 for the lattice widened on 29 Oct. Total about 42 min.

| # | step | takes | check, and what is normal |
|---|---|---|---|
| 1 | `dark30`: 30 s dark at 0 pA | 0.5 min, plus Jin Yong's prompt | `Overview dark30`: resting Vm. Healthy on 28 Oct: −58 mV. The degraded rerun read −41 mV. Then Jin Yong **types the holding current** (`Holding current recorded as "..."`). The current monitor reads a constant 2–4 pA offset; that is normal. |
| 2 | `p1`: Protocol 1 coarse map | 6.2 min | `Overview p1`: Vm (−67 mV on 28 Oct). `FLAG p1: Vm drifted ±D mV` fires above 5 mV of drift. |
| 3 | `p1_peak` | seconds | `DECISION p1_peak/aim: frame F at col/row [c r]`. Uses the 6-px grid peak when ≥ 2 mV, else the 12-px grid (FLAG). Below 2 mV on both grids the session **stops** (`FLAG session stopped at p1_peak: Protocol 1 peak ... below the floor ...`): no usable response, so Jin Yong decides (new cell, or a `'peak_frame'` override). |
| 4 | `p2a`: Protocol 2 stage 1, the sweeps | 8.3 min | `Overview p2a`. Watch for the drift FLAG. |
| 5 | `readout1` | seconds | `DECISION readout1/placement_pd: placement PD X deg (rule mean ...)`. Normal: rule `mean`, 28 and 56 °/s within a few degrees (1.1° on 28 Oct). Rule `28_only` means they differ by more than 10°: note it. `usable_pd: stage 2 plays the 22 bars` is normal. `FLAG stage 2 fallback: all 88 bars as today (N rows, M min)` makes step 6 last M min instead of 2: tell Jin Yong. |
| 6 | `p2b`: stage 2, 22 bars × 3 | 2.2 min | |
| 7 | `readout` | seconds | `Read-out (...): ... center col/row [c r], ... center peak X mV`, then `DECISION readout/center`. Normal: `read-out M6 center; tables shifted by [0 0] px`. `FLAG DECISION readout/center: Protocol 1 peak ...` means the center failed a sanity check (more than 6 px from the Protocol 1 peak, more than 4 px off the bar set, or a center peak below 3 mV). The session continues aimed at the Protocol 1 peak: tell Jin Yong and note it. |
| 8 | `build` | seconds | `Prot3b: 26 frames, 26 blocks, 7.6 min; Protocol 4: ~450 frames, ~1780 flashes, ~15 min.` `FLAG Prot3b build:` or `FLAG Protocol 4 build:` list squares off the screen. A few is fine; most is not. |
| 9 | `p3b`: Prot3b, ±0.5 pA via AO2 | 7.8 min | `Analysis p3b: R_gray G GOhm; pairs kept k/8; dR ...`. 28 Oct: 13 GΩ, 6/8, ON center −5.7 GΩ. `FLAG p3b: only k of 8 pairs pass the baseline QC` means fewer than 4 pairs kept: the cell drifted. The cell is not dropped for drift; drift is a result. |
| 10 | `p4`: Protocol 4, combined map | ~15.4 min | `Analysis p4: lattice ON peak X mV; ... passes kept k/6`. 28 Oct: 17.4 mV, 6/6. `FLAG p4: k of 6 passes within 5 mV of pass 1` means the cell drifted. |
| – | end | | `Arena off (allOff): panels dark.` then `Session complete.` A `FLAG arena not turned off` means Jin Yong must switch the arena off by hand. |

After every step, report: the step, its duration against the table, Vm, and any DECISION or FLAG in plain words.
Example: *"p2a done, 8.3 min, Vm −70 mV, steady. Placement PD 271°, both speeds agree. Stage 2 next (2 min)."*

**Signs the cell is going** (run 2 on 28 Oct is the example): Vm climbing step after step (−41 → −25 mV), drift
FLAGs on consecutive steps, R_gray falling (13 → 5.4 GΩ), Prot3b pairs and Protocol 4 passes rejected. Say so
early and give the numbers. Whether to continue is Jin Yong's call.

## 6. When something goes wrong

| what you see | what it means | tell Jin Yong |
|---|---|---|
| `Step X FAILED: <msg>` | the session stopped; the arena is cleaned up | Quote the message and the stack from `Auto_T4_v1_*.log`. **The same command resumes** (`Auto_T4_v1('<TAG>', 'L')`): run_session re-runs the failed step and continues. If the same step fails twice, or the message points at the code (an index or frame error, like `stage-2 frames outside 2..89` on 28 Oct), do not resume blindly: record it as a bug (section 7). |
| `FLAG session stopped at <step>: ...` (e.g. after `p1_peak`) | run_session stopped on purpose (the step is `interrupted`), e.g. no usable Protocol 1 response | New cell, or a `'peak_frame'` override, at Jin Yong's choice. |
| `run_session:controller`, `startDisplay did not confirm`, Host errors | rig communication problem | Stop. Do not suggest workarounds on the hardware. Record it. Resuming is Jin Yong's decision. |
| Ctrl-C in their MATLAB | the step is interrupted, with controller cleanup | The same command resumes. |
| the cell is lost | | `Auto_T4_v1('<TAG>', 'L', 'abandon', true)` closes the session. The next cell gets a **new tag**. |
| `FLAG arena not turned off` | | Switch the arena off by hand. |
| no new log line for much longer than the step's duration | possible hang | Check `session_state.json` (still `running`?) and the tail of `Auto_T4_v1_*.log`. Tell Jin Yong; do not act. |

## 7. Your record: `C:\Users\labadmin\session_results\<TAG>\sidecar_notes.md`

Keep it up to date as you go. Include: the tag, side, date; the holding current; one line per step (start time,
duration, Vm, DECISIONs and FLAGs); everything Jin Yong tells you about the cell (seal, access, noise, anything
unusual at the rig); every failure, with the exact message, what Jin Yong did and the result. Anything that looks
like a code bug goes under a `## Bugs` heading: step, exact message, log excerpt, the stack from
`Auto_T4_v1_*.log`, and what happened when the session was resumed.

## 8. After `Session complete.`

Build the results page (MATLAB is allowed now; it only reads the session folders):
```powershell
$o = 'C:\Users\labadmin\session_results\<TAG>'
$cmd = "addpath(genpath('C:\matlabroot\GitHub\nested_RF_stimulus_session_ui\src')); session_results_page({'D:\nested_RF_sessions\<TAG>'}, '$o\page', 'Title', '<TAG>')"
Start-Process 'C:\Program Files\MATLAB\R2023a\bin\matlab.exe' -ArgumentList '-batch',"`"$cmd`"" -NoNewWindow -Wait `
  -RedirectStandardOutput "$o\page_stdout.txt" -RedirectStandardError "$o\page_stderr.txt"
```
`$o\page\index.html` is one self-contained file. It holds the numbers table, the standard figure of each protocol,
the Protocol 4 detail figure and the step overviews. To compare with earlier cells, pass several folders:
`{'D:\nested_RF_sessions\JY1028a', 'D:\nested_RF_sessions\<TAG>'}`. Add `'Labels', {...}` for short names.

Then give Jin Yong a short summary in this order: recording quality (Vm per step, drift), the Protocol 1 peak, the
placement PD, the center and its distance to the Protocol 1 peak, Prot3b (R_gray, pairs kept, the three ΔR
read-outs), Protocol 4 (lattice peak, passes kept, the trailing hyperpolarization). Compare with the
28 Oct healthy cell, and list every FLAG and failure. Finish `sidecar_notes.md` with the same summary.

## 9. Background, if Jin Yong asks why

- The sequence, its parameters and the a priori read-outs: `C:\matlabroot\ephys-standard-analysis\docs\briefs\JY_PROTOCOL_HANDOVER_20261026.md`.
  The data behind each choice: `PROTOCOL_REDESIGN_20261025.md` in the same folder.
- Changes since that hand-over:
  - **Prot3b pair QC:** the two blocks of the *same* current in each ABBA quartet must be within 2 mV of each
    other, and every block within 5 mV of the first block of its current (agreed 29 Oct). The hand-over's
    wording would reject every pair.
  - **Protocol 4 lattice window:** now s −14…+18 px (it was −12…+16; widened 29 Oct, about 1.1 min more).
- Conventions: 1 px = 1.25°; s is positive toward the PD (the trailing side); v is 90° clockwise from the PD.
