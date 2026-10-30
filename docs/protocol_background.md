# Protocol background

Design history and rationale for the nested receptive field (RF) and direction selectivity
(DS) protocols. For instructions on running them, see the
[README](../README.md) and the
[Quarto documentation](https://leburnett.github.io/reiser-documentation/Ephys/ephys_nested_rf.html).

## Origin

The nested RF protocol was designed by LE Burnett in 2024 for patch electrophysiology
experiments run by JY Park in the Reiser Lab at HHMI's Janelia Research Campus. It was based
on the protocol used in Gruntman et al. 2019.

## Why two protocols

The two protocols split a problem that cannot be solved in one pass: mapping a receptive
field at high spatial resolution requires knowing roughly where the receptive field is
first, and probing direction selectivity requires knowing the cell's preferred contrast.

**Protocol 1** determines the rough receptive field location and the preferred contrast of
the recorded neuron. It is presented from a set of pre-made `.g4p` files, organised by which
side of the G4 arena the stimulus appears on. After presentation, processing scripts analyse
the data and generate plots showing the response of the cell to dark and bright flashes at
different positions on the arena. From these plots the user reads off the `peak_frame` — the
frame number that elicited the greatest response.

**Protocol 2** takes that `peak_frame` as input and probes direction selectivity using moving
bar stimuli, while measuring receptive field structure at higher spatial resolution using
small flashing squares centred on the identified location.

## Protocol 2 stimulus design

Generated de novo on every run by `generate_protocol2()`, except the bar sweep and bar flash
stimuli, which use pre-made bar patterns stored in `results/patterns/protocol2`. To change
the parameters of the bar sweep or bar flash stimuli, those pre-made patterns must be
modified directly — regenerating the protocol will not change them.

| Component | Design |
|---|---|
| 4 px flash grid | 196 flashes, 14 x 14 grid, 50% overlap, 30 px crop around centre |
| 6 px flash grid | 100 flashes, 10 x 10 grid, 50% overlap, 33 px crop around centre |
| Bar patterns | Cropped and centred on the identified `[x, y]` coordinate, 30 px crop |
| Flash timing | 160 ms flash, 440 ms inter-flash interval (600 ms per flash) |
| Pixel intensities | background 4, off 0, on 15 |

## Reference

Gruntman, Romani & Reiser (2019).
