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
addpath('C:\matlabroot\GitHub\nested_RF_stimulus_reaim\src\protocol_generation'); Auto_T4_v2('<TAG>', 'L')      % with the re-aim
addpath('C:\matlabroot\GitHub\nested_RF_stimulus_session_ui\src\protocol_generation'); Auto_T4_v1('<TAG>', 'L')  % fallback
```
Resume a session folder with the version that started it; a fallback cell gets a new tag.

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
   - The ON and OFF squares whose center lies within `reaim_window_px` of the prediction.
   - Protocol 1's timing: 200 ms flash, 150 ms gray. Each repeat is one row: 1 s gray, then the squares in a seeded
     random order.
   - `reaim_reps` repeats, with 2 s of gray before and after.

   On the 6-px grid this is 16–25 squares × ON/OFF × 2, about 35 s. `rapid_p1_peak` then finds the ON peak in the
   window, and that peak replaces `out.p1_peak` for everything after: the Protocol 2 aim, the read-out's `P1Peak`
   and the center check. The old peak is kept in `out.p1_peak.before_reaim`.
4. **Decision `p1_confirm/confirm`.**
   - **GOOD:** the found row is within `reaim_tolerance_px` of the predicted row, the peak is not on the window's
     outer row or column, it has at least the margin on both sides, and it moved the predicted way.
   - **Anything else:** `CHECK`, flagged with the reasons. The found peak is still used, because it was measured
     after the move.
   - **No ON response of at least `p1_peak_floor_mV` in the window:** the session stops. The message gives the
     stage reading to go back to and suggests `'reaim_reset', true`, which makes the step ask again on resume.
5. **The figure and the choice** (decision `p1_confirm/gate`). The confirmation map is drawn on the arena and saved
   as `plots\p1_confirm_reaim_<time>.png`. It shows the window's ON ΔV, the margins, the peak before the move, and the
   predicted and found positions. The operator then chooses:
   - **Enter**: continue to Protocol 2.
   - **c**: repeat the check.
   - **a**: re-aim again from the found peak, starting from the current stage reading. Offered only when it is
     still near the edge.
   - **s**: stop. Resuming runs the check again.

   After `c` or `a`, the main loop sets those steps back to pending and runs them again; the status table shows
   this.

The whole re-aim adds about 1.5–2 min to an edge cell: the rotation, then 40 s of check, then the choice.

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
| `reaim_window_px` | [12 12] | ±2 squares of the 6-px grid around the prediction, enough for an error in the angle or in the degrees-per-pixel of up to about 9 px (one-and-a-half squares) and still well inside the screen. 16–25 positions, about 6 % of the 6-px grid. |
| `reaim_reps` | 2 | Protocol 1 uses 4 repeats. Two give a mean map that finds a response of several mV, at about 35 s. |
| `reaim_flash_s`, `reaim_gray_s` | 0.2, 0.15 | Protocol 1's own timing, so the ΔV is comparable with the Protocol 1 map. |
| `reaim_tolerance_px` | 6 | One 6-px square: the found peak lies on the grid square nearest the true RF. |
| `reaim` | `'auto'` | `'off'` logs and flags the edge but asks nothing. |
| `reaim_reset` | false | After a stop at `p1_confirm` (or to redo a re-aim): `reaim` and `p1_confirm` become pending again and the Protocol 1 peak is restored. Rotate the stage back first. |
| `operator_answers` | {} | Scripted answers to the prompts, for the tests. |

## Rig checklist (first test on the real arena)

**Code**

1. `addpath('C:\matlabroot\GitHub\nested_RF_stimulus_reaim\src\protocol_generation'); Auto_T4_v2('<TAG>', 'L')`
   (`Auto_T4_v2` checks for `feat/arena-reaim` by default).
2. Preflight first: the same command with `'check_only', true`. Passed on 29 Oct for both v2 and v1.
3. Fallback, if v2 misbehaves:
   `addpath('C:\matlabroot\GitHub\nested_RF_stimulus_session_ui\src\protocol_generation'); Auto_T4_v1('<NEW TAG>', 'L')`.

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

8. At the prompt, rotate by the suggested angle and type the angle actually moved. Expect `DECISION
   p1_confirm/confirm: GOOD: found col/row ..., predicted ...` within about 6 px.
9. If the flags say `moved the OPPOSITE way ... reaim_sign is wrong`, or the session stops at `p1_confirm` with no
   response:
   - Rotate back to the starting reading.
   - Resume with `'reaim_reset', true, 'reaim_sign', <the other sign>`.
   - Change the default in `run_session_v2.m` once it is confirmed.
10. Compare the measured shift with the prediction: rows moved (in `confirm`'s inputs, found minus `before_xy`)
    against `actual / 1.25`. A consistent ratio away from 1 means `reaim_deg_per_px` should change. Note it in the
    rig log after 2–3 cells.
11. After the session, check `p2a`'s line `aimed at col/row [c r]`: it must be the confirmed position.

**Not covered yet**

- The confirmation's ΔV is not compared with the Protocol 1 peak automatically. Both are in the inputs.
- The arena viewer shows `p1_confirm` like the other stimulus steps (one case added).
- No rotation of the base back to 0 is asked at the end of the session. Pass `stage_angle_deg` for the next cell,
  or turn the base back by hand.
