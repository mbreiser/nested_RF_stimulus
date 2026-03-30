function [canvas, shift_amount] = shift_traces_to_canvas(traces_11xN, centroid_idx)
% SHIFT_TRACES_TO_CANVAS  Place 11-position traces onto a 15-position canvas.
%
%   [CANVAS, SHIFT_AMOUNT] = SHIFT_TRACES_TO_CANVAS(TRACES_11XN, CENTROID_IDX)
%   places an 11×N matrix of bar flash traces onto a 15×N NaN-padded canvas,
%   shifted so that the response-weighted centroid aligns to the canvas center
%   (column 8). This is a linear shift — no circular wrapping.
%
%   The nearest-neighbor shift is computed as:
%       shift_amount = 6 - round(centroid_idx)
%   where centroid_idx is the fractional RF centroid position (1–11, center=6).
%
%   Canvas layout (15 columns):
%       Column 8 = aligned center (0°)
%       Unshifted cell: positions 1–11 → canvas columns 3–13
%       Cell with shift=+2: positions 1–11 → canvas columns 5–15
%       Cell with shift=-1: positions 1–11 → canvas columns 2–12
%
%   Degree labels for 15 columns (2.5° spacing):
%       −17.5°, −15.0°, −12.5°, ..., 0° (col 8), ..., +15.0°, +17.5°
%
%   INPUTS:
%     traces_11xN - 11×N matrix of baseline-subtracted bar flash traces
%                   (11 positions × N timepoints), ordered ND to PD
%     centroid_idx - Fractional centroid position (1–11), where 6 = center
%
%   OUTPUTS:
%     canvas       - 15×N matrix with NaN padding where no data exists
%     shift_amount - Integer shift applied (positive = shift toward PD side)
%
%   EXAMPLE:
%     % Cell with centroid at position 4.2 (shifted +2 toward PD)
%     [canvas, shift] = shift_traces_to_canvas(data, 4.2);
%     % shift = +2, data placed in canvas rows 5–15
%
%   See also COMPUTE_RF_CENTROID, ANALYZE_SINGLE_EXPERIMENT_MR

    n_canvas = 15;
    n_orig   = 11;
    n_pad    = (n_canvas - n_orig) / 2;  % = 2

    N = size(traces_11xN, 2);  % number of timepoints

    % Nearest-neighbor shift
    shift_amount = 6 - round(centroid_idx);

    % Place the 11 original positions onto the canvas
    % Unshifted: start at column n_pad+1 = 3
    start_row = n_pad + 1 + shift_amount;

    % Initialise canvas with NaN
    canvas = NaN(n_canvas, N);

    % Place data (clamping to valid range should not be needed for ±2 shifts,
    % but guard against unexpected centroid values)
    row_start = max(1, start_row);
    row_end   = min(n_canvas, start_row + n_orig - 1);

    % Corresponding source indices
    src_start = row_start - start_row + 1;
    src_end   = src_start + (row_end - row_start);

    canvas(row_start:row_end, :) = traces_11xN(src_start:src_end, :);

end
