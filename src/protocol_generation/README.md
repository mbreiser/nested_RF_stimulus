# src/protocol_generation

Builds a protocol 2 experiment centred on the receptive-field location found in protocol 1,
and runs it on the G4 arena.

## Entry points

| Run this | What it does | Calls |
|---|---|---|
| `generate_protocol2()` | The whole flow: prompt, generate, assemble, run | `get_input_parameters`, `create_protocol2`, `run_protocol2` |
| `find_inv_peak_frame.m` | Inverse peak frame for the opposite-contrast experiment | — |
| `export_to_google_sheets.m` | Pushes experiment metadata to the recording log | — |

## Inputs

Prompted for in a dialog, not passed as arguments:

| Input | Meaning |
|---|---|
| `peak_frame` | Frame number that gave the largest response in protocol 1 |
| Arena side | Which side of the G4 arena protocol 1 was presented on |
| Fly age, strain | Recorded into `currentExp.mat` as metadata |

`peak_frame` is converted to screen coordinates `[x, y]`, which also determines the ON/OFF
contrast preference used for the flash stimuli.

## Outputs

A timestamped folder `yyyy_MM_dd_HH_mm` under
`<matlabroot>\G4_Protocols\nested_RF_protocol2\`, containing:

| Item | Contents |
|---|---|
| `Patterns/` | 4 px flash grid (196 flashes, 14 x 14), 6 px grid (100 flashes, 10 x 10), cropped bar patterns centred on `[x, y]` |
| `Functions/` | Position functions for the flash and bar stimuli |
| `currentExp.mat` | Pattern and function ordering, plus the metadata above |
| `Log/` | Recorded data, written during the run |

## Prerequisites

- Protocol 1 has been run and a `peak_frame` identified.
- `G4_Display_Tools` is installed and configured.
- The arena is connected and calibrated.

## Gotchas

- **The bar sweep and bar flash stimuli are not generated here.** They use pre-made patterns
  in `results/patterns/protocol2`; to change their parameters, edit those patterns directly.
  Everything else is regenerated on every run.
- Stimulus timing and intensities are fixed in the code: 160 ms flash, 440 ms inter-flash
  interval, pixel intensities background 4 / off 0 / on 15. See
  [docs/protocol_background.md](../../docs/protocol_background.md).
