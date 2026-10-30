# scripts

Top-level scripts you run directly: the manuscript figures and the batch results table.

## Entry points

| Run this | What it does | Calls |
|---|---|---|
| `generate_manuscript_fig_main.m` | Main manuscript figure | `generate_manuscript_fig('main')` |
| `generate_manuscript_fig_supp.m` | Supplementary figure | `generate_manuscript_fig('supp')` |
| `build_batch_results.m` | Builds the pooled results table across cells | — |

Each figure wrapper is a few lines: it sets `DATA_ROOT` and calls the combiner.

## Before you run

Set `DATA_ROOT` at the top of the wrapper script to your local copy of the `ttl_1DRF` tree,
and `CIRCSTAT_PATH` in `build_batch_results.m` to your CircStat installation.

## Full detail

[MANUSCRIPT_FIGURES.md](MANUSCRIPT_FIGURES.md) documents the figure pipeline: which function
draws which panel, the expected data layout under `DATA_ROOT`, the arguments of each
sub-figure function, and what was done by hand in Illustrator.

## Gotchas

- Outputs are timestamped, so repeated runs accumulate rather than overwrite.
- `build_batch_results.m` needs the Circular Statistics Toolbox on the MATLAB path.
