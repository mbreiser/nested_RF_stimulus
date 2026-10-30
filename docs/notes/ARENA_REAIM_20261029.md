# Arena re-aim after Protocol 1 (29 Oct 2026)

Branch `feat/arena-reaim` (cut from `feat/session-progress-plots` at `ddfac2b`), checkout
`C:\matlabroot\GitHub\nested_RF_stimulus_reaim`.
- **New version:** launcher `src/protocol_generation/Auto_T4_v2.m` and sequencer `run_session_v2.m`, a copy of
  `run_session.m` with the steps `reaim` and `p1_confirm`.
- **Fallback:** `Auto_T4_v1.m` and `run_session.m` are unchanged (byte-identical to `feat/session-progress-plots`),
  and stay the last working version.
- **Tests:** `src/test_protocols/test_run_session_reaim.m` and `test_run_session_v2_dry.m`. `test_run_session_dry.m`
  still tests the unchanged `run_session.m`.
- **Sidecar's view:** `docs/EXPERIMENT_SIDECAR.md`, sections 0 and 5a.

```matlab
% with the re-aim
addpath('C:\matlabroot\GitHub\nested_RF_stimulus_reaim\src\protocol_generation')
Auto_T4_v2('<TAG>', 'L')
% fallback (restart MATLAB first if in doubt)
addpath('C:\matlabroot\GitHub\nested_RF_stimulus_session_ui\src\protocol_generation')
Auto_T4_v1('<TAG>', 'L')
```
Resume a session folder with the version that started it; a fallback cell gets a new tag. Type both lines for
every call: since 30 Oct `Auto_T4_v2` takes every `_reaim` folder off the path when it returns, fails or is stopped,
and `Auto_T4_v1.m` is not on this branch (see "30 Oct changes" below).

## 30 Oct changes (after JY1030a, 29 Oct)

**What happened.** The Protocol 1 ON peak was 9.8 mV at SNR 10 (frame 161, col/row [37.5 45.5]). The stage moved
−30°. The ±12-px confirmation window (2 repeats, ON and OFF) returned 3.0 mV at SNR 1.7, a noise checkerboard of
−2.3 to +3 mV, and was judged GOOD. Protocol 2 was aimed at noise, the bar responses were poor, and the session was
stopped. The cause: the noise per square of a 2-repeat mean is about 1.4 mV (single flash about 2 mV), so the
largest of 25 squares of noise alone is about 3 mV, and the check had only the 2 mV floor. Separately, at 17:12
Jin Yong could not start the v1 fallback: the `_reaim` path left over from the v2 run made MATLAB run this branch's
copy of `Auto_T4_v1.m`, which refuses to start outside `feat/session-progress-plots`.

**Changes** (M. Reiser's spec of 30 Oct):
1. The confirmation scans the whole Protocol 1 area: every ON square of the Protocol 1 pattern, no OFF (6-px LHS:
   128 squares, columns 17–112), 2 repeats, Protocol 1's timing, 95.6 s. The peak comes from the full map. The
   prediction stays an overlay, and its offset is reported as a geometry flag, not a failure. The window is kept
   as `'reaim_confirm', 'window'`.
2. The QC: (a) SNR ≥ `reaim_min_snr` (4), as `rapid_p1_peak`'s `snr_mean`; (b) peak ≥ `reaim_min_frac` (0.5) × the
   Protocol 1 ON peak before the move, and ≥ `p1_peak_floor_mV`; (c) not on the scan's outer row or column. (a)
   and (b) pass and (c) passes: GOOD. (a) and (b) pass, (c) fails: CHECK, still used. (a) or (b) fails: NO
   RESPONSE. Protocol 2 is not aimed there, and the choices are c (repeat), r (rotate back, redo the re-aim), p
   (rotate back, go on with the peak from before the move), u (use the scan's peak anyway; added later on 30 Oct at
   M. Reiser's request, so that a clear but weaker peak never costs the cell; flagged `operator_accept`), s (stop).
   The SNR and the ratio are in `decisions.json`, in the figure's title and in the prompt. The prediction never
   decides anything: a response is accepted wherever the scan finds it.
3. The path: `Auto_T4_v1.m` is removed from this branch, so the fallback can only run from `_session_ui`.
   `Auto_T4_v2` restores the path at every exit (return, error, Ctrl-C) without any `_reaim` folder, the typed
   `addpath` included. Checked on 30 Oct: after `Auto_T4_v2(..., 'check_only', true)` and after an error, no `_reaim`
   folder is on the path. With both checkouts on the path (the `_reaim` one first), `which Auto_T4_v1` names the
   `_session_ui` copy, and its `check_only` passes.

Not changed, pending M. Reiser's bench check of 30 Oct: `reaim_sign`, `reaim_deg_per_px`, and any column shift
(there is no column-shift model).

**Why.** When the Protocol 1 peak is too close to the top or bottom edge, the later protocols lose much of their
stimulus off the screen. Jin Yong used to move the display by hand, or the cell was lost. The arena sits on a
swiveling base with a range of ±20°. The sequencer now detects the edge case and asks for a rotation. It then
confirms the new position with a short Protocol 1, so these cells can be recorded. Cells away from the edge are not
affected, and nothing is asked.

## What happens

1. **`reaim`**, right after `p1_peak`. The room is the screen above and below the peak's center: `row − 0.5`
   above and `48.5 − row` below, where the row is the median pixel of the peak square. If both are at least
   `reaim_margin_px`, the step logs `edge_check` and does nothing else.
2. If either side has less room, the operator sees a boxed prompt with this suggested rotation:
   `suggested = reaim_sign × (reaim_target_row − row) × reaim_deg_per_px`. It is limited to ±`reaim_limit_deg` on
   the stage scale (starting from `stage_angle_deg`) and rounded toward 0 in steps of `reaim_step_deg`. The
   operator rotates the base and presses Enter (or types `s` to skip), then types the angle actually moved. The
   predicted new row is `row + reaim_sign × actual / reaim_deg_per_px`; the column is unchanged.
3. **`p1_confirm`** runs only after a rotation. It plays a package in `<session>\p1_confirm\`:
   - Protocol 1's own square pattern, copied as pattern 1 (the `.mat` and the `.pat`). The frame numbers are the
     same as in Protocol 1, so `rapid_p1_peak` reads it unchanged. The Protocol 1 package itself is never edited.
   - Every ON square of that pattern (`'reaim_confirm', 'full'`, the default), or, with `'window'`, the ON and OFF
     squares whose center lies within `reaim_window_px` of the prediction.
   - Protocol 1's timing: 200 ms flash, 150 ms gray. Each repeat is one row: 1 s gray, then the squares in a seeded
     random order.
   - `reaim_reps` repeats, with 2 s of gray before and after.

   On the 6-px LHS grid the full scan is 128 squares × 2 repeats: 95.6 s. (The window is 16–25 squares × ON/OFF × 2,
   about 35 s.) `rapid_p1_peak` then finds the ON peak of the map.
4. **Decision `p1_confirm/confirm`**, the QC (30 Oct):
   - **A response** needs (a) SNR ≥ `reaim_min_snr` and (b) peak ≥ max(`p1_peak_floor_mV`, `reaim_min_frac` × the
     Protocol 1 ON peak before the first move). It then replaces `out.p1_peak` for everything after: the Protocol 2
     aim, the read-out's `P1Peak` and the center check. The old peak is kept in `out.p1_peak.before_reaim`.
   - **GOOD:** a response that is (c) not on the outer row or column of the scan. **CHECK:** a response on it
     (flagged; still used).
   - **NO RESPONSE:** (a) or (b) fails, or the read-out's peak is not a square of the scan. `out.p1_peak` is not
     changed.
   - **Geometry flags** (never a failure): an offset from the prediction above `reaim_tolerance_px` in col or row,
     a move the opposite way to the prediction, still near the edge.
5. **The figure and the choice** (decision `p1_confirm/gate`). The scan's ON map is drawn on the arena and saved as
   `plots\p1_confirm_reaim_<time>.png`, with the scanned area, the margins, the peak before the move, and the
   predicted and found positions. The title gives the mV, the SNR and the ratio. The operator then chooses:
   - after GOOD or CHECK: **Enter** continue to Protocol 2; **c** repeat the scan; **a** re-aim again from the found
     peak with the current stage reading (offered only when it is still near the edge); **s** stop (resuming runs the
     scan again).
   - after NO RESPONSE: **c** repeat the scan; **r** rotate the stage back to the reading before the move (then
     Enter), and the re-aim is asked again; **p** rotate back (then Enter) and go on with the peak from before the
     move; **u** use the scan's peak anyway (the operator judges it real; Protocol 2 aimed there, flagged
     `operator_accept`, `out.p1_peak.accepted_by_operator`); **s** stop (the message gives the stage reading and
     `'reaim_reset', true`). Enter is not a choice here.

   After `c`, `a` or `r`, the main loop sets those steps back to pending and runs them again; the status table shows
   this.

The whole re-aim adds about 2.5–3 min to an edge cell: the rotation, then 1.6 min of scan, then the choice.

In `run_session_v2`, a step that stops after logging a decision returns `S.stop_request` to the main loop, instead
of raising the stop error inside the step. Raising it inside the step lost the decision from the state file, and
`decisions.json` was then overwritten on resume. `run_session.m` keeps the old behavior (its `p1_peak` stops).

## Parameters and why

| option | default | why |
|---|---|---|
| `reaim_margin_px` | 14 | Protocol 2 crops its bars and flashes to 30 px around the center (`px_crop_bar` = `px_crop_flash` = 30 in `generate_protocol2_package`, i.e. 15 px each side; it already warns outside rows 16..33). The Prot3b core sites reach 14 px: the trailing site at s +12 plus half of a 4-px square (`place_sites` in `rapid_readout.m`). The ±14 "outside" sites are mirrored when off screen, so they do not set the margin. 14 keeps the 6-px grid's rows 15.5 and 33.5, where Protocol 2 loses at most one edge row, and re-aims rows 9.5 / 39.5 and beyond. On the 12-px grid it re-aims rows 6.5 and 42.5. Any value from 9.5 to 15 gives the same classification on the 6-px grid. |
| `reaim_target_row` | 24.5 | The middle of rows 1–48: 24 px of room on each side. |
| `reaim_deg_per_px` | 1.25 | The arena's pixel pitch. The cylinder makes it approximate away from the equator; the confirmation measures the real shift. |
| `reaim_sign` | +1 | +1 means a positive rotation on the stage scale moves the RF **down** the screen (toward row 48). **Unconfirmed**: find it with the checklist below. |
| `reaim_limit_deg` | 20 | The stage range, ±20°. |
| `stage_angle_deg` | 0 | The stage reading when the session starts. The clamp is applied to the absolute reading. |
| `reaim_step_deg` | 0.5 | A readable step on the scale. Rounding toward 0 never goes past the range. |
| `reaim_confirm` | `'full'` | Since 30 Oct, every ON square of the Protocol 1 pattern: the peak is found wherever the rotation put it, so a wrong angle, sign or column shift cannot hide it. ON only, because the aim uses the ON peak; this halves the time. `'window'` is the 29 Oct check. |
| `reaim_window_px` | [12 12] | `'window'` only: ±2 squares of the 6-px grid around the prediction. 16–25 positions. |
| `reaim_reps` | 2 | Protocol 1 uses 4 repeats. Two give 95.6 s for the full scan. With 2 repeats the noise per square is about 1.4 mV (single flash about 2 mV, JY1030a), which is why the QC below is needed. |
| `reaim_min_snr` | 4 | SNR of the ON peak as `rapid_p1_peak`'s `snr_mean`: peak / (noise / √repeats), the noise being the RMS over the ON squares of the SD across repeats. For noise alone, the largest of 128 squares is about 2.6 SD (the noise-only test: 3.6 mV, SNR 2.6; JY1030a: SNR 1.7). A peak at 4 SD from noise alone has a probability of about 0.4 % over 128 squares. With JY1030a's noise, an RF at the full Protocol 1 amplitude (9.8 mV) gives an SNR of about 7 with 2 repeats. One at half of it (4.9 mV) is at the limit of both tests (SNR about 3.5). |
| `reaim_min_frac` | 0.5 | The RF should respond about as well after the move as before. Half the Protocol 1 ON peak catches a small but clean response (high SNR at low noise) that is not the RF. The 2 mV floor (`p1_peak_floor_mV`) still applies. |
| `reaim_flash_s`, `reaim_gray_s` | 0.2, 0.15 | Protocol 1's own timing, so the ΔV is comparable with the Protocol 1 map. |
| `reaim_tolerance_px` | 6 | One 6-px square. A larger offset (col or row) from the prediction is a geometry flag: the measured peak is used either way. |
| `reaim` | `'auto'` | `'off'` logs and flags the edge but asks nothing. |
| `reaim_reset` | false | After a stop at `p1_confirm` (or to redo a re-aim): `reaim` and `p1_confirm` become pending again and the Protocol 1 peak is restored. Rotate the stage back first. |
| `operator_answers` | {} | Scripted answers to the prompts, for the tests. |

## Rig checklist (first test on the real arena)

**Code**

1. Two lines (`Auto_T4_v2` checks for `feat/arena-reaim` by default):
   ```matlab
   addpath('C:\matlabroot\GitHub\nested_RF_stimulus_reaim\src\protocol_generation')
   Auto_T4_v2('<TAG>', 'L')
   ```
2. Preflight first: the same two lines with `'check_only', true`. Passed on 30 Oct for both v2 and v1.
3. Fallback, if v2 misbehaves (restart MATLAB first if in doubt; `which Auto_T4_v1` must name `_session_ui`):
   ```matlab
   addpath('C:\matlabroot\GitHub\nested_RF_stimulus_session_ui\src\protocol_generation')
   Auto_T4_v1('<NEW TAG>', 'L')
   ```

**See the prompt once without a cell (2 min)**

4. Run a dry run with a top-edge peak. Frame 171 of the default 6-px LHS pattern is the ON square at rows 7–12,
   cols 47–52 (checked):
   `run_session_v2(fullfile(tempdir, 'reaim_try'), 'dry_run', true, 'auto_accept', false, 'p1_package', '<the 6pxONLY LHS package>', 'peak_frame', 171, 'steps', {'p1', 'p1_peak', 'reaim', 'p1_confirm'})`.
   It asks for Enter and the angle typed, as on the rig.

**Find the rotation sign (no cell, 5 min, arena on)**

5. **Rows.** Play Protocol 1 frame 2 of the 6-px pattern. It is the dark square at rows 1–6, cols 17–22. Check that
   it is at the **top** of the screen as seen from the fly's position. If it is at the bottom, the row convention is
   upside down: reverse the result of the next step.
6. **Direction.** Turn the base to +5° on its scale, and watch which way the screen in front of the fly moves:
   - **Screen moves up**: a point fixed in the fly's view now falls on a lower part of the screen, so the RF moves
     toward row 48. Keep `reaim_sign` +1.
   - **Screen moves down**: pass `'reaim_sign', -1` (to `Auto_T4_v2`, which passes it to `run_session_v2`).

   Turn the base back to 0.
7. Note which way positive is on the scale in the rig log. Each cell starts at a known reading: 0, or pass
   `'stage_angle_deg', X`.

**First edge cell**

8. At the prompt, rotate by the suggested angle and type the angle actually moved. Expect, after 1.6 min,
   `DECISION p1_confirm/confirm: GOOD: found col/row ... (X mV, SNR Y (needs 4), Z x the peak before the move
   (needs 0.5)), predicted ...`. Geometry flags (offset, opposite way) do not stop anything: note them.
9. **NO RESPONSE** (Protocol 2 is not aimed anywhere yet): look at the figure and at the Vm in `Overview p1_confirm`.
   - A clear peak, only weaker than before (or just under the SNR line): `u` uses it.
   - The cell is going: `s`, or `p` to record the rest at the old angle.
   - The RF probably left the screen (wrong sign): `r`, rotate back when asked, and at the new prompt turn the other
     way and type the angle moved with its sign. Or `s`, then resume with `'reaim_reset', true, 'reaim_sign', <the
     other sign>`.
   - Change the default `reaim_sign` in `run_session_v2.m` only after the bench check.
10. Compare the measured shift with the prediction: rows moved (in `confirm`'s inputs, `found_xy` minus `before_xy`)
    against `actual / 1.25`, and any column shift (`offset_px(1)`). A consistent ratio away from 1 means
    `reaim_deg_per_px` should change. Note it in the rig log after 2–3 cells.
11. After the session, check `p2a`'s line `aimed at col/row [c r]`: it must be the confirmed position (or the old
    one after `p`).

**Not covered yet**

- No model of a column shift with the rotation; the full scan finds the peak wherever it is.
- The arena viewer shows `p1_confirm` like the other stimulus steps (one case added).
- No rotation of the base back to 0 is asked at the end of the session. Pass `stage_angle_deg` for the next cell,
  or turn the base back by hand.
