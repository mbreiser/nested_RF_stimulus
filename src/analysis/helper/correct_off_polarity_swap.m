function [bar_data, was_corrected] = correct_off_polarity_swap( ...
    bar_data, pattern_order, func_order, date_str, is_off, func_pair)
% CORRECT_OFF_POLARITY_SWAP  Fix direction swap from dark-bar patterns 9/10.
%
%   [BAR_DATA, WAS_CORRECTED] = CORRECT_OFF_POLARITY_SWAP(BAR_DATA, ...
%       PATTERN_ORDER, FUNC_ORDER, DATE_STR, IS_OFF)
%   corrects for a stimulus polarity error in the pre-bar-flash (summer 2025)
%   experiments. Pattern files 0009 and 0010 used dark bars (bar_0) instead
%   of bright bars (bar_15), which caused func=3 and func=4 response
%   directions to appear swapped for OFF cells.
%
%   [BAR_DATA, WAS_CORRECTED] = CORRECT_OFF_POLARITY_SWAP(..., FUNC_PAIR)
%   specifies which function pair to swap. Default is [3 4] (slow/28 dps).
%   Use [5 6] for medium/56 dps bars.
%
%   The correction swaps bar_data rows for the specified func pair of
%   experiment patterns 9 and 10, so that the LUT direction assignments
%   become correct.
%
%   APPLICABILITY:
%     - Only applied to OFF cells (IS_OFF == true)
%     - Only applied to experiments before Oct 15, 2025 (the date the
%       pattern files were corrected on the arena SD card)
%     - ON cells are unaffected because bright/dark polarity does not
%       change the direction of their excitatory response
%     - Post-Oct-15 OFF cells (e.g. 2025_10_22_13_10) already have
%       corrected patterns on the arena and do NOT need the swap
%
%   INPUTS:
%     bar_data      - Nx(R+1) cell array from parse_bar_data_pre_bf. Rows
%                     correspond to bar directions in presentation order
%                     (matching metadata pattern/func ordering)
%     pattern_order - Full 1xM pattern order from currentExp.mat
%     func_order    - Full 1xM function order from currentExp.mat
%     date_str      - Experiment date string, 'YYYY_MM_DD_HH_MM'
%     is_off        - true if this is an OFF/T5 cell
%     func_pair     - (Optional) [f_fwd f_rev] function numbers to swap.
%                     Default: [3 4] (slow/28 dps). Use [5 6] for 56 dps.
%
%   OUTPUTS:
%     bar_data      - Corrected cell array (unchanged if no swap needed)
%     was_corrected - true if the correction was applied
%
%   BACKGROUND (commit 82badc7c):
%     Commit "update 15 and 16 patterns - correct direction" (Oct 15 2025)
%     fixed patterns 0015/0016 from bar_0 -> bar_15, but patterns 0009-0014
%     remain with bar_0 in the repository. The arena SD card patterns were
%     updated around the same time. Empirical verification on all 23
%     experiments confirms: all 7 pre-Oct-15 OFF cells need the swap, and
%     all ON cells + the 1 post-Oct-15 OFF cell do not.
%
%   See also PARSE_BAR_DATA_PRE_BF, VERIFY_LUT_DIRECTIONS,
%            BATCH_ANALYZE_PRE_BAR_FLASH

    was_corrected = false;

    if nargin < 6 || isempty(func_pair), func_pair = [3 4]; end

    % --- Check whether correction is needed ---
    if ~is_off
        return;
    end

    if ~needs_correction(date_str)
        return;
    end

    % --- Build bar index for the target func pair ---
    f_fwd = func_pair(1);
    f_rev = func_pair(2);
    target_mask  = (func_order == f_fwd | func_order == f_rev);
    target_pats  = pattern_order(target_mask);
    target_funcs = func_order(target_mask);

    % --- Swap forward <-> reverse for patterns 9 and 10 ---
    for pat = [9, 10]
        row_f3 = find(target_pats == pat & target_funcs == f_fwd, 1);
        row_f4 = find(target_pats == pat & target_funcs == f_rev, 1);

        if isempty(row_f3) || isempty(row_f4)
            warning('correct_off_polarity_swap:missingRows', ...
                'Cannot find both func %d/%d for pattern %d. Skipping.', ...
                f_fwd, f_rev, pat);
            continue;
        end

        % Swap all columns (per-cycle traces + mean)
        temp = bar_data(row_f3, :);
        bar_data(row_f3, :) = bar_data(row_f4, :);
        bar_data(row_f4, :) = temp;
    end

    was_corrected = true;

end


%% ======================================================================

function tf = needs_correction(date_str)
% NEEDS_CORRECTION  Check if experiment date is before Oct 15, 2025.

    parts = sscanf(date_str, '%d_%d_%d_%*d_%*d');

    if numel(parts) < 3
        warning('correct_off_polarity_swap:dateParseError', ...
            'Could not parse date from "%s". Assuming no correction.', date_str);
        tf = false;
        return;
    end

    exp_date = datenum(parts(1), parts(2), parts(3));
    cutoff   = datenum(2025, 10, 15);
    tf       = exp_date < cutoff;

end
