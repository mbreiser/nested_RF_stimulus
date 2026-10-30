# Experiment sidecar: supervising an Auto_T4_v2 (or Auto_T4_v1) session

For a Claude session running on the rig PC (REISERLAB-EPHYS) next to Jin Yong's experiment. Written 29 Oct 2026
from the first real sessions (JY1028a/b, 28 Oct). **Claude: read all of this before you do anything.**

## 0. Which version is running

There are two launchers. Ask Jin Yong which one, or read the first line of the `Auto_T4_v*_<time>.log` in the session folder.

| | **Auto_T4_v2** (from 29 Oct: with the arena re-aim) | **Auto_T4_v1** (the fallback, the last working version) |
|---|---|---|
| checkout | `C:\matlabroot\GitHub\nested_RF_stimulus_reaim` (branch `feat/arena-reaim`) | `C:\matlabroot\GitHub\nested_RF_stimulus_session_ui` (branch `feat/session-progress-plots`) |
| command | `addpath('C:\matlabroot\GitHub\nested_RF_stimulus_reaim\src\protocol_generation'); Auto_T4_v2('<TAG>', 'L')` | `addpath('C:\matlabroot\GitHub\nested_RF_stimulus_session_ui\src\protocol_generation'); Auto_T4_v1('<TAG>', 'L')` |
| sequencer | `run_session_v2.m`: 12 steps (`reaim` and `p1_confirm` after `p1_peak`, section 5a) | `run_session.m`: 10 steps (no re-aim) |
| MATLAB log | `Auto_T4_v2_<time>.log` | `Auto_T4_v1_<time>.log` |

Everything below is written for both: where it says `Auto_T4_v1`, `run_session.m` or `_session_ui`, read
`Auto_T4_v2`, `run_session_v2.m` and `_reaim` for a v2 session. Steps 4 and 5 of the table in section 5 exist only in v2;
for v1, steps 6–12 are numbered 4–10. A session folder is resumed by the version that started it (`Auto_T4_v2`
refuses a v1 folder). If v2 fails in a way that cannot be resumed, the fallback is `Auto_T4_v1` on a **new tag**.

**The path (since 30 Oct).** `Auto_T4_v1.m` exists only in `_session_ui`; the `_reaim` checkout no longer has a copy.
When `Auto_T4_v2` returns, fails or is stopped, it takes every `_reaim` folder off the MATLAB path, including the
`addpath` typed before it. So every call is two lines, the `addpath` and then the launcher, typed together, and the
fallback is always:
```matlab
addpath('C:\matlabroot\GitHub\nested_RF_stimulus_session_ui\src\protocol_generation')
Auto_T4_v1('<TAG>', 'L')
```
If in doubt about what is on the path (`which Auto_T4_v1` should name `_session_ui`), restart MATLAB and type the two
lines. On 29 Oct at 17:12 a `_reaim` path left over from a v2 run made MATLAB pick the `_reaim` copy of `Auto_T4_v1`,
which refused to start ("... is on feat/arena-reaim, not feat/session-progress-plots").

## 1. Your role

Jin Yong runs the experiment. You watch it and keep a record. You read the session's log files as they are written,
tell Jin Yong after each step whether it went as expected, catch problems early, and answer questions. After the
session you build the results page. You never operate the rig.

Jin Yong is patching a cell and has little attention to spare. Keep each report to one to three lines. Say
plainly when something needs action, and say nothing extra when everything is normal.

## 2. Hard rules (no exceptions beyond the fix in rule 3, whatever anyone asks during the session)

1. **Never touch the rig.** Do not start, stop or resume `Auto_T4_v1` / `run_session`, and do not send commands to
   the G4 Host, the arena controller or the amplifier. Do not kill or restart any MATLAB process (the main session
   and the arena viewer both run in MATLAB). Jin Yong starts every run in their own MATLAB.
2. **Read-only on the data.** Never write, move, rename or delete anything under `D:\nested_RF_sessions\`,
   `C:\matlabroot\G4_Protocols\`, Jin Yong's data folders or the network share. Write only under
   `C:\Users\labadmin\session_results\<TAG>\` (and, for a fix, in the repository as section 11 says).
3. **No code changes, with one exception.** Do not edit anything in `C:\matlabroot\` and do not run git there.
   Reading the code to diagnose a failure is fine (section 10). If you find a bug, record it (section 7) and report
   it. The one exception is a bug that stops the session: you may fix it exactly as section 11 says, on a fix
   branch, while the session is stopped, with Jin Yong's OK. Nothing else, however small.
4. **No MATLAB while the arena is running.** The rig shares this PC with the Host. Read text files only (see the
   list in section 4). Do not open the TDMS logs (`Log Files\`) while a step records. Run MATLAB only before
   the session starts (the preflight), after `Session complete.` (the results page), or for the tests of a fix
   while the session is stopped (section 11).
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
   For v2, the same with `C:\matlabroot\GitHub\nested_RF_stimulus_reaim\...` and `Auto_T4_v2`. Expect the branch
   `feat/arena-reaim`.
4. Jin Yong then runs, in their own MATLAB, after break-in in current clamp at 0 pA:
   ```matlab
   addpath('C:\matlabroot\GitHub\nested_RF_stimulus_reaim\src\protocol_generation'); Auto_T4_v2('<TAG>', 'L')
   ```
   The fallback (the last working version, no re-aim), two lines; restart MATLAB first if in doubt:
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
- `---- k/12  <step>  (<label>, about N min) ----`: a step starts (12 steps since the re-aim was added; `reaim` and
  `p1_confirm` are steps 4 and 5)
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
| 4 | `reaim`: arena re-aim check | seconds; 1–2 min when it fires | Normal for most cells: `DECISION reaim/edge_check: no re-aim: A px above and B px below the peak (>= 14)`, and step 5 is skipped (`Step p1_confirm skipped: no re-aim`). **When it fires** (`DECISION reaim/edge_check: RE-AIM: peak D px from the top/bottom edge (< 14)`), the session **waits for Jin Yong** at a boxed `RE-AIM THE ARENA` prompt in their MATLAB (see section 5a). That wait is expected, not a hang. |
| 5 | `p1_confirm`: Protocol 1 ON scan of the whole Protocol 1 area at the new angle, then the figure and Jin Yong's choice | about 1.6 min plus the choice, only after a rotation | `DECISION p1_confirm/confirm: GOOD: found col/row [c r] (X mV, SNR Y (needs 4), Z x the peak before the move (needs 0.5)), predicted [c r]; Protocol 2 is aimed there` is a good confirmation. `FLAG ... CHECK` (peak on the scan's outer row or column) and `FLAG ... NO RESPONSE` (Protocol 2 **not** aimed there; Jin Yong must choose) are in section 5a. |
| 6 | `p2a`: Protocol 2 stage 1, the sweeps | 8.3 min | `Overview p2a`. Watch for the drift FLAG. After a re-aim, the `Protocol 2 package generated for peak frame F ... aimed at col/row [c r]` line must show the **confirmed** position, not the Protocol 1 one. |
| 7 | `readout1` | seconds | `DECISION readout1/placement_pd: placement PD X deg (rule mean ...)`. Normal: rule `mean`, 28 and 56 °/s within a few degrees (1.1° on 28 Oct). Rule `28_only` means they differ by more than 10°: note it. `usable_pd: stage 2 plays the 22 bars` is normal. `FLAG stage 2 fallback: all 88 bars as today (N rows, M min)` makes step 8 last M min instead of 2: tell Jin Yong. |
| 8 | `p2b`: stage 2, 22 bars × 3 | 2.2 min | |
| 9 | `readout` | seconds | `Read-out (...): ... center col/row [c r], ... center peak X mV`, then `DECISION readout/center`. Normal: `read-out M6 center; tables shifted by [0 0] px`. `FLAG DECISION readout/center: Protocol 1 peak ...` means the center failed a sanity check (more than 6 px from the Protocol 1 peak, more than 4 px off the bar set, or a center peak below 3 mV). The session continues aimed at the Protocol 1 peak: tell Jin Yong and note it. |
| 10 | `build` | seconds | `Prot3b: 26 frames, 26 blocks, 7.6 min; Protocol 4: ~450 frames, ~1780 flashes, ~15 min.` `FLAG Prot3b build:` or `FLAG Protocol 4 build:` list squares off the screen. A few is fine; most is not. |
| 11 | `p3b`: Prot3b, ±0.5 pA via AO2 | 7.8 min | `Analysis p3b: R_gray G GOhm; pairs kept k/8; dR ...`. 28 Oct: 13 GΩ, 6/8, ON center −5.7 GΩ. `FLAG p3b: only k of 8 pairs pass the baseline QC` means fewer than 4 pairs kept: the cell drifted. The cell is not dropped for drift; drift is a result. |
| 12 | `p4`: Protocol 4, combined map | ~15.4 min | `Analysis p4: lattice ON peak X mV; ... passes kept k/6`. 28 Oct: 17.4 mV, 6/6. `FLAG p4: k of 6 passes within 5 mV of pass 1` means the cell drifted. |
| – | end | | `Arena off (allOff): panels dark.` then `Session complete.` A `FLAG arena not turned off` means Jin Yong must switch the arena off by hand. |

After every step, report: the step, its duration against the table, Vm, and any DECISION or FLAG in plain words.
Example: *"p2a done, 8.3 min, Vm −70 mV, steady. Placement PD 271°, both speeds agree. Stage 2 next (2 min)."*

### 5a. The arena re-aim (steps 4 and 5; Auto_T4_v2 only, added 29 Oct)

**Why.** Protocol 2 crops its bars and flashes to 30 px around the Protocol 1 peak, and the Prot3b sites reach
14 px from it. A peak within 14 px of the top or bottom edge (on the 6-px grid: rows 3.5, 9.5, 39.5, 45.5) would
leave much of those stimuli off the screen. Before 29 Oct, Jin Yong moved the display by hand or gave up the cell.
Now the sequencer asks for a rotation of the swiveling base instead, so these cells can be recorded. Cells whose peak
is not near an edge are not affected, and nothing is asked.

**The prompt** (Jin Yong's MATLAB, boxed `RE-AIM THE ARENA`). It gives the peak's col/row, its distance to the edge,
and an angle such as `Rotate the stage by +18.5 deg (from +0 to +18.5 deg on its scale)`. It also gives the expected
move: `the RF moves 14.8 px DOWN (toward row 48), to row 24.3`. The angle re-centers the peak on row 24.5:
(24.5 − row) × 1.25° × `reaim_sign`, rounded to 0.5° and limited to the stage's ±20°. When the limit cuts it,
the prompt says `CLAMPED`. Jin Yong rotates the base and presses Enter, then types the angle actually moved
(Enter accepts the suggestion). Typing `s` skips the re-aim, and the session continues at the edge peak. A long
wait at this prompt is normal. Do not report it as a hang.

**What is recorded** (`decisions.json`):
- `reaim/edge_check`: the room above and below the peak, and the margin (14 px).
- `reaim/rotation`: the angle wanted, suggested and actually moved, whether it was clamped, the stage reading
  before and after, the sign convention, and the predicted new col/row. Flags: `moved X deg, suggested Y`,
  `clamped ...`, `operator_skip`, `stage_at_limit`.
- `p1_confirm/confirm`: the scan's ON peak (frame, col/row, mV), its SNR and noise, the ratio to the Protocol 1 peak
  before the move, which QC tests passed (`qc`), the predicted col/row and the offset.
- `p1_confirm/gate`: Jin Yong's choice after the figure (below).

**The check** (since 30 Oct) is an ON-only scan of the **whole Protocol 1 area** at the new angle: every ON square
of the Protocol 1 pattern (6-px LHS grid: 128 squares, columns 17–112, rows 1–48), 2 repeats at 200 + 150 ms, 1 s
of gray before each repeat and 2 s before and after. It takes **95.6 s (1.6 min)**. The peak comes from the full
map, so it is found even when the predicted angle was wrong. The 29 Oct version, a ±12-px window of ON and OFF
squares around the prediction (about 35 s), is still available as `'reaim_confirm', 'window'`.

**The QC: the peak must be real.** On 29 Oct (JY1030a) the window check reported `GOOD` at 3.0 mV with SNR 1.7.
That was noise: the mean of 2 repeats has about 1.4 mV of noise, so the largest of many squares of noise alone is
about 3 mV, over the 2 mV floor. Protocol 2 was then aimed at noise. Now the scan's ON peak counts as a response only
when:
- (a) its SNR (computed as `rapid_p1_peak` computes `snr_mean`) is at least `reaim_min_snr` (4), and
- (b) it is at least `reaim_min_frac` (0.5) × the Protocol 1 ON peak before the move (and at least 2 mV).

A response is **GOOD** when (c) it is not on the scan's outer row or column. It is **CHECK** when it is: it is still
used, but the RF may lie beyond the screen or the grid. If (a) or (b) fails, the verdict is **NO RESPONSE** and
Protocol 2 is **not** aimed there. The offset from the predicted position is only a geometry flag, never a failure
(the angle per pixel and any column shift are being measured on the bench, 30 Oct).

**The figure** opens and is saved as `plots\p1_confirm_reaim_<time>.png`. It shows the arena with the scan's ON ΔV
map, the scanned area (dashed), the dotted 14-px margins, the Protocol 1 peak before the move (gray square), the
predicted position (red ×) and the found peak: green ○ GOOD, orange CHECK, red "largest ON (not a response)" for
NO RESPONSE. The title and the prompt give the mV, the SNR and the ratio. Jin Yong then chooses in the command window.

After GOOD or CHECK:
- **Enter**: continue, with Protocol 2 aimed at the found peak.
- **c**: repeat the scan, stage unmoved.
- **a**: re-aim again from the found peak, with the stage where it is. Offered only when the found peak is still
  near the edge.
- **s**: stop the session. The same command resumes it: the scan runs again, then the same choice.

After NO RESPONSE (Enter is not a choice here):
- **c**: repeat the scan, stage unmoved.
- **r**: rotate the stage back to the reading before the move (the prompt names it, then Enter), and the re-aim is
  asked again.
- **p**: rotate the stage back to that reading (then Enter) and go on with the Protocol 1 peak from before the move.
- **s**: stop. To be asked the re-aim again on resume: rotate back and add `'reaim_reset', true` (and
  `'reaim_sign', -1` if the base turned the wrong way). The same command without it repeats the scan.

`Running again from step ...` in the log follows `c`, `a` or `r`. The status table then shows the step pending
again. Waiting at this choice is normal.

The confirmation package and its event ledger are in `<TAG>\p1_confirm\`, and the read-out is in
`<TAG>\readout\p1_confirm\`. The Protocol 1 peak from before the move is kept in `session_state.json` under
`out.p1_peak.before_reaim`.

**A good confirmation:** `GOOD`, SNR well over 4, a ratio near 1. Report it in one line, for example:
*"Re-aimed +18.5°, RF found at row 21.5 (predicted 24.3), 7.4 mV, SNR 11, 0.9 × before; Protocol 2 aimed there."*

**Flags to report** (in brackets after the decision):
- `geometry: found [dc dr] px (col, row) from the prediction ...` or `geometry: moved the OPPOSITE way ...: reaim_sign is
  wrong`. The aim is still right, because it uses the measured peak. Note the numbers for the bench calibration.
- `peak on the outer row or column of the scan`: the verdict is CHECK. The aim is uncertain, so watch the
  `readout/center` check later.
- `still near the edge`: the rotation was not enough, for example because it was clamped.
- **NO RESPONSE** with `SNR x < 4` and/or `peak x mV < y mV (r x the ... before the move; needs 0.5 x)`. Either the
  cell is going (check the Vm in `Overview p1_confirm`), or the RF went off the screen. Tell Jin Yong the SNR and the
  ratio. Usually the best choice is **p** (the cell is still recorded at the old angle, with the edge peak), or
  **r** to try the rotation again.

**Signs the cell is going** (run 2 on 28 Oct is the example): Vm climbing step after step (−41 → −25 mV), drift
FLAGs on consecutive steps, R_gray falling (13 → 5.4 GΩ), Prot3b pairs and Protocol 4 passes rejected. Say so
early and give the numbers. Whether to continue is Jin Yong's call.

## 6. When something goes wrong

| what you see | what it means | tell Jin Yong |
|---|---|---|
| `Step X FAILED: <msg>` | the session stopped; the arena is cleaned up | Quote the message and the stack from `Auto_T4_v1_*.log`. **The same command resumes** (`Auto_T4_v1('<TAG>', 'L')`): run_session re-runs the failed step and continues. If the same step fails twice, or the message points at the code (an index or frame error, like `stage-2 frames outside 2..89` on 28 Oct), do not resume blindly: find the cause in the code (section 10), record it as a bug (section 7), and offer Jin Yong a fix (section 11). |
| `FLAG session stopped at <step>: ...` (e.g. after `p1_peak`) | run_session stopped on purpose (the step is `interrupted`), e.g. no usable Protocol 1 response | New cell, or a `'peak_frame'` override, at Jin Yong's choice. |
| `FLAG DECISION p1_confirm/confirm: NO RESPONSE ...`, then the choice c / r / p / s | the re-aim scan has no real response (SNR or amplitude too low): the RF left the screen, or the cell is going | Section 5a. Give the SNR and the ratio. `p` goes on at the old angle with the old peak; `r` tries the rotation again. |
| `FLAG session stopped at p1_confirm: no response in the re-aim scan; stopped by the operator ...` | Jin Yong chose `s` after NO RESPONSE | Rotate the base back, then resume with `'reaim_reset', true` (and `'reaim_sign', -1` if the base turned the wrong way). |
| `Unrecognized function or variable 'Auto_T4_v1'`, or `Auto_T4_v1: ... is on feat/arena-reaim` | the `_session_ui` folder is not first on the path | Section 0: type the `_session_ui` `addpath` line, then `Auto_T4_v1`. If in doubt, restart MATLAB first. |
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
`Auto_T4_v1_*.log`, the cause you found in the code (file, function, line), and what happened when the session
was resumed. For a fix (section 11), also: the diff, Jin Yong's two OKs in their words, the test results, the
commit and branch, and which steps ran on the fixed code.

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
  - **Arena re-aim** (29 Oct, Auto_T4_v2, section 5a). Parameters and the rig test:
    `docs\notes\ARENA_REAIM_20261029.md` in the `_reaim` checkout.
- Conventions: 1 px = 1.25°; s is positive toward the PD (the trailing side); v is 90° clockwise from the PD.

## 10. Where the code is

The rig runs this checkout, `C:\matlabroot\GitHub\nested_RF_stimulus_session_ui` (branch
`feat/session-progress-plots`, or the day's fix branch, section 11). `Auto_T4_v1` puts its `src\` first on the path
and checks that the files below come from it. The stack trace in `Auto_T4_v1_*.log` names the file and line; start
there.

| what | where (under `src\`) |
|---|---|
| the command Jin Yong runs: checks, the viewer, then run_session | `protocol_generation\Auto_T4_v1.m` (in `_session_ui` only) |
| the sequencer: state, resume, decisions, the log, the step loop | `protocol_generation\run_session.m` (main function) |
| `dark30`, `p1` | `run_session.m`: `step_dark30`, `step_p1`, `rows_from_g4p`, `run_rows` |
| `p1_peak` | `run_session.m`: `step_p1_peak`, `p1_square_patterns`, `frame_to_xy`; the read-out `analysis\protocol2\rapid_p1_peak.m` |
| `p2a` (stage 1) | `run_session.m`: `step_p2_stage1`, `generate_protocol2_package`; `protocol_generation\create_protocol2.m`, `run_protocol2.m`; `stimulus_generation\posn_fn\generate_bar_pos_fns.m`, `generate_bar_flash_pos_fns.m` |
| `readout1`, `readout` | `run_session.m`: `step_readout1`, `step_readout`, `call_readout`, `read_readout`; `analysis\protocol2\rapid_readout.m` (the Protocol 4 window: `lattice_window`), `compute_m6_centroid.m`, `pipeline\plot_rapid_readout.m`, `pipeline\plot_rapid_stage1.m`, `pipeline\fast_tdms_folder2struct.m` |
| `p2b` (stage 2) | `run_session.m`: `step_p2_stage2`, `build_stage2_package` (the 28 Oct `frames outside 2..89` bug was here) |
| `build` | `run_session.m`: `step_build`, `build_prot3b_package`, `build_prot4_package`, `pats_from_frame_table`, `pfn_sections`, `write_pattern`, `write_function` |
| `p3b`, `p4` | `run_session.m`: `step_run_package`, `run_rows`, `controller_cleanup` |
| end of session | `run_session.m`: `arena_off` |
| overviews and the analysis lines (they never stop the session) | `run_session.m`: `step_overview`; `protocol_generation\session_step_overview.m`; `analysis\session_plots\analyze_prot3b.m`, `analyze_p4.m` |
| results page, offline analysis | `analysis\session_plots\session_results_page.m`, `analyze_session.m`, `plot_p4_detail.m`, `session_figure.m`, `load_step_events.m`, `step_snips.m` |
| progress window, arena viewer | `protocol_generation\SessionProgress.m`, `arena_viewer.m` |
| mock controller for the dry runs | `protocol_generation\SessionMockController.m` |
| tests (the ones section 11 runs) | `test_protocols\test_run_session_dry.m`, `test_run_session_replay_jy.m`, `test_run_session_p1_aim.m`, `test_session_plots.m` (`test_run_session_joint.m` only runs on a Mac) |
| **v2 only**: the launcher and the sequencer | `protocol_generation\Auto_T4_v2.m`, `run_session_v2.m` (a copy of `run_session.m` with the re-aim; the function names above are the same in it) |
| **v2 only**: `reaim`, `p1_confirm` | `run_session_v2.m`: `step_reaim`, `step_p1_confirm` (the QC and the gate after the figure), `confirm_snr`, `build_p1_confirm_package`, `plot_reaim_confirm`, `reaim_reset`; `Auto_T4_v2.m`: `restore_path` (takes the `_reaim` folders off the path at every exit) |
| **v2 only**: tests | `test_run_session_reaim.m`, `test_run_session_v2_dry.m` (run them with the four above, which test v1's `run_session.m`) |

Outside this repository (read only, never edit, even in a fix):
- `C:\matlabroot\G4_Display_Tools\`: arena control (`PControl_Matlab\controller\PanelsController.m`) and the
  Host configuration. A trace that ends there is a rig problem (section 6), not a code fix.
- `C:\matlabroot\G4_Protocols\`: the Protocol 1 packages and the stage-1 data.
- `C:\matlabroot\GitHub\ephys-standard-analysis\`: the hand-over and redesign briefs (`docs\briefs\`), the Python
  read-out (`tools\rapid_prot2_locations.py`), and the reference tables the dry test reads (`tools\out\`).
- `C:\matlabroot\GitHub\nested_RF_stimulus\` and the other checkouts: not what the rig runs.

## 11. Fixing a bug that stops the session

Only for a bug that **stops the session**: `Step X FAILED` with a stack trace into this repository's `src\`, where
resuming would fail again. Everything else waits until after the session and goes under `## Bugs`: a FLAG, a failed
overview or analysis line (the session goes on; the figures can be redone offline), a wrong-looking number, a
trace into `G4_Display_Tools` or the Host.

What a fix may and may not do:
- It makes the code do what the hand-over specifies, with the smallest change that does it: an index, a bounds
  check, a missing field, a wrong file name.
- It does not change any stimulus parameter (sizes, durations, positions, the lattice window, currents, repeats),
  any threshold or any decision rule. If the only fix would, it is not yours to make: tell Jin Yong, who decides
  whether to end the session.
- It touches only files under `src\` in this repository. Prefer `run_session.m` and the read-outs over
  `run_protocol2.m`, `create_protocol2.m` and `generate_protocol2.m`, which Jin Yong also uses by hand.
- Budget about 15 minutes from the failure to the resume. The cell does not wait longer. If you are not there by
  then, say so and stop; Jin Yong decides.

The procedure. Do the steps in order and do not skip any:
1. **Confirm the session is stopped.** Jin Yong's MATLAB prompt is back (`Auto_T4_v1` has returned) and
   `session_state.json` shows the step `failed`, with no step `running`. Never edit while a step runs: MATLAB picks
   up a changed file the next time it calls it, even in the same run.
2. **Diagnose.** Find the line from the stack trace and explain the cause in two or three sentences.
3. **Propose.** Show Jin Yong the exact change (the lines before and after) and say which steps will run on the
   fixed code. Ask: *"Fix it on the fix branch and run the tests?"* Go on only on a clear yes, and write their
   words into the notes.
4. **Branch and edit.** Use the git that ships with GitHub Desktop. The checkout must have no changes before you
   start. Git's CRLF warnings can make it return exit code 255 even when it worked, so read its output.
   ```powershell
   $g = 'C:\Users\labadmin\AppData\Local\GitHubDesktop\app-3.5.3\resources\app\git\cmd\git.exe'
   $r = 'C:\matlabroot\GitHub\nested_RF_stimulus_session_ui'
   & $g -C $r status --short --branch           # nothing but the branch line; note the branch
   & $g -C $r switch -c fix/sidecar-<yyyyMMdd>  # once per day; if it exists and is checked out, stay on it
   ```
   Make the edit, then show `& $g -C $r diff`.
5. **Test.** All four tests together take about 80 s. They write only under the temp folder and use the mock
   controller, never the arena:
   ```powershell
   $o = 'C:\Users\labadmin\session_results\<TAG>\fix_tests'; New-Item -ItemType Directory -Force $o | Out-Null
   $t = 'C:\matlabroot\GitHub\nested_RF_stimulus_session_ui\src\test_protocols'
   $p = foreach ($n in 'test_run_session_dry', 'test_run_session_replay_jy', 'test_run_session_p1_aim', 'test_session_plots') {
     Start-Process 'C:\Program Files\MATLAB\R2023a\bin\matlab.exe' -ArgumentList '-batch',"`"run('$t\$n.m')`"" -NoNewWindow -PassThru `
       -RedirectStandardOutput "$o\$n.txt" -RedirectStandardError "$o\$n.err" }
   $p | Wait-Process
   Select-String -Path "$o\*.txt" -Pattern 'check\(s\) failed' | ForEach-Object { "$($_.Filename): $($_.Line)" }
   Get-ChildItem "$o\*.err" | ForEach-Object { "$($_.Name): $((Get-Content $_.FullName) -join ' | ')" }
   ```
   Pass: each of the four reports `0 check(s) failed.`, and the `.err` files are empty, except
   `test_run_session_dry.err`, which normally holds one line: `run_session: controller cleanup after an
   interruption (10 steps, 10 ok)`. If anything else fails, the fix does not go in. Commit it anyway, so it can be
   reviewed, with `NOT TESTED OK` in the message. Switch back to the branch you noted in step 4, and tell Jin Yong
   that the session can only resume on the old code (which will probably fail again) or end.
6. **Commit, never push.** Commit only the files you changed:
   ```powershell
   & $g -C $r commit -m "Sidecar fix during <TAG>: <what>" -m "<step> failed with <message>; cause <...>. Tests: 4/4 pass. Approved by Jin Yong at <HH:mm>." -- <changed files>
   & $g -C $r log --oneline -1
   ```
   Do not push, merge, rebase or reset, and do not touch the other branches or worktrees.
7. **Hand back.** Ask Jin Yong: *"Tests pass; resume on the fix?"* On yes, they resume in their own MATLAB, with the
   fix branch named, because `Auto_T4_v1` checks the branch:
   ```matlab
   Auto_T4_v1('<TAG>', 'L', 'branch', 'fix/sidecar-<yyyyMMdd>')
   ```
   Watch the step that failed, closely. Record the commit, and that every step from the failed one on ran on it.

For the rest of the day, every `Auto_T4_v1` call, including the preflight of the next cell, needs
`'branch', 'fix/sidecar-<yyyyMMdd>'`. A second fix that day goes on the same branch. Never switch back to
`feat/session-progress-plots` yourself, except to undo a fix that failed its tests (step 5). The fix is reviewed
and merged after the experiment, by the lab, not by you. If a later day's preflight stops with
`... is on fix/sidecar-..., not feat/session-progress-plots`, an earlier fix has not been reviewed yet. Tell Jin Yong,
show them `& $g -C $r log --oneline feat/session-progress-plots..HEAD`, and let them decide whether to run on that
branch.
