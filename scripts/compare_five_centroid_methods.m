% COMPARE_FIVE_CENTROID_METHODS  Compare 5 RF center-finding methods on all 25 cells.
%
%   Loads the saved amplitude profiles (peak_amplitudes) from centering
%   summary .mat files and applies 5 different center-finding methods to
%   each cell's 1x11 amplitude vector. No pipeline re-run needed.
%
%   Methods:
%     M1  Assumed center (always position 6)
%     M2  Robust peak (position of maximum amplitude)
%     M3  Flanked peak (peak validated by top-4 neighbors on both sides)
%     M4  Full depolarizing centroid (center-of-mass over all A>0 positions)
%     M5  FWHM bump centroid (center-of-mass restricted to half-max bump)
%
%   Output: comparison table with discrepancy highlighting.

%% Setup
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));

data_root   = '/Users/reiserm/Documents/ttl_1DRF';
results_dir = fullfile(data_root, 'population_results');

%% Load saved centering summaries
S_on  = load(fullfile(results_dir, 'on_cells_centering_summary.mat'));
S_off = load(fullfile(results_dir, 'off_cells_centering_summary.mat'));

all_cells = [S_on.on_cells;  S_off.off_cells];
all_infos = [S_on.cell_infos; S_off.cell_infos];
n_total   = size(all_cells, 1);

%% Apply 5 methods to each cell
%  Store results in arrays
m1 = nan(n_total, 1);  % assumed center
m2 = nan(n_total, 1);  % robust peak
m3 = nan(n_total, 1);  % flanked peak (NaN = no bump)
m4 = nan(n_total, 1);  % full depolarizing centroid
m5 = nan(n_total, 1);  % FWHM bump centroid
m3_flag = false(n_total, 1);  % true = no bump

for i = 1:n_total
    ci = all_infos{i};
    A  = max(ci.peak_amplitudes, 0);  % threshold negatives to 0

    % M1: Assumed center
    m1(i) = 6.00;

    % M2: Robust peak
    m2(i) = method_robust_peak(A);

    % M3: Flanked peak
    [m3(i), m3_flag(i)] = method_flanked_peak(A);

    % M4: Full depolarizing centroid
    m4(i) = method_full_centroid(A);

    % M5: FWHM bump centroid
    m5(i) = method_fwhm_centroid(A);
end

%% Convert to degrees offset from center
deg = @(idx) (idx - 6) * 2.5;

%% Print comparison table
fprintf('\n');
fprintf('=============================================================================================================================================\n');
fprintf('  5-METHOD RF CENTER COMPARISON  (all values are position indices, 1-11; position 6 = 0 deg)\n');
fprintf('=============================================================================================================================================\n');
fprintf('%-26s %-12s   M1     M2     M3      M4      M5    MaxD   Flag\n', ...
    'Folder', 'Group');
fprintf('%-26s %-12s  (=6)  (peak) (flnk)  (full)  (FWHM)\n', '', '');
fprintf('%s\n', repmat('-', 1, 105));

n_flagged = 0;
flagged_idx = [];

for i = 1:n_total
    % MaxDelta: max pairwise difference among M2-M5
    vals = [m2(i), m4(i), m5(i)];
    if ~m3_flag(i)
        vals = [vals, m3(i)];
    end
    max_delta = max(vals) - min(vals);

    % Flag marker
    flag_str = '';
    if max_delta > 1.0
        flag_str = '  <<<';
        n_flagged = n_flagged + 1;
        flagged_idx(end+1) = i; %#ok<SAGROW>
    end

    % M3 display
    if m3_flag(i)
        m3_str = '  NB  ';
        flag_str = '  <<<';  % always flag no-bump cells
        if max_delta <= 1.0
            n_flagged = n_flagged + 1;
            flagged_idx(end+1) = i; %#ok<SAGROW>
        end
    else
        m3_str = sprintf('%5d ', m3(i));
    end

    fprintf('%-26s %-12s  %4.0f   %4d   %s %6.2f  %6.2f  %5.2f%s\n', ...
        all_cells{i, 1}, all_cells{i, 2}, ...
        m1(i), m2(i), m3_str, m4(i), m5(i), max_delta, flag_str);
end

fprintf('%s\n', repmat('-', 1, 105));
fprintf('%d of %d cells flagged (MaxDelta > 1.0 or M3 = no bump)\n\n', ...
    numel(unique(flagged_idx)), n_total);

%% Print integer shift comparison
fprintf('=== INTEGER SHIFT COMPARISON (shift = 6 - round(centroid)) ===\n');
fprintf('%-26s %-12s  S1  S2  S3  S4  S5   Agree?\n', 'Folder', 'Group');
fprintf('%s\n', repmat('-', 1, 80));

n_disagree = 0;
for i = 1:n_total
    s1 = 6 - round(m1(i));
    s2 = 6 - round(m2(i));
    s4 = 6 - round(m4(i));
    s5 = 6 - round(m5(i));
    if m3_flag(i)
        s3_str = 'NB';
        shifts = [s1 s2 s4 s5];
    else
        s3 = 6 - round(m3(i));
        s3_str = sprintf('%+d', s3);
        shifts = [s1 s2 s3 s4 s5];
    end
    all_agree = all(shifts == shifts(1));
    if ~all_agree
        n_disagree = n_disagree + 1;
        agree_str = ' ***';
    else
        agree_str = '';
    end

    fprintf('%-26s %-12s  %+d  %+d  %s  %+d  %+d%s\n', ...
        all_cells{i,1}, all_cells{i,2}, ...
        s1, s2, s3_str, s4, s5, agree_str);
end
fprintf('\n%d of %d cells have disagreeing integer shifts (marked ***)\n\n', ...
    n_disagree, n_total);

%% Amplitude profiles for flagged cells
if ~isempty(flagged_idx)
    flagged_idx = unique(flagged_idx);
    fprintf('=== AMPLITUDE PROFILES FOR FLAGGED CELLS ===\n');
    for k = 1:numel(flagged_idx)
        i = flagged_idx(k);
        ci = all_infos{i};
        A = max(ci.peak_amplitudes, 0);

        fprintf('\n  %s (%s)\n', all_cells{i,1}, all_cells{i,2});
        fprintf('  Pos:  ');
        for p = 1:11, fprintf('%6d', p); end
        fprintf('\n');
        fprintf('  Amp:  ');
        for p = 1:11, fprintf('%6.1f', A(p)); end
        fprintf('\n');

        % Show rank ordering
        [~, rank_order] = sort(A, 'descend');
        ranks = zeros(1, 11);
        for r = 1:11, ranks(rank_order(r)) = r; end
        fprintf('  Rank: ');
        for p = 1:11
            if A(p) == 0
                fprintf('     -');
            else
                fprintf('%6d', ranks(p));
            end
        end
        fprintf('\n');
        fprintf('  M2=pos%d  M3=%s  M4=%.2f  M5=%.2f\n', ...
            m2(i), ...
            ternary(m3_flag(i), 'NB', sprintf('pos%d', m3(i))), ...
            m4(i), m5(i));
    end
    fprintf('\n');
end

%% Summary statistics by group
groups = {'on_control', 'on_ttl', 'off_control', 'off_ttl'};
fprintf('=== SUMMARY BY GROUP ===\n');
fprintf('%-12s   n    M4-M5 (mean+/-SD)   M2-M5 (mean+/-SD)   M3 no-bump\n', 'Group');
fprintf('%s\n', repmat('-', 1, 75));

for g = 1:numel(groups)
    mask = strcmp(all_cells(:,2), groups{g});
    n_g  = sum(mask);

    diff_45 = m4(mask) - m5(mask);
    diff_25 = m2(mask) - m5(mask);
    n_nb    = sum(m3_flag(mask));

    fprintf('%-12s  %2d    %+5.2f +/- %.2f       %+5.2f +/- %.2f        %d/%d\n', ...
        groups{g}, n_g, ...
        mean(diff_45), std(diff_45), ...
        mean(diff_25), std(diff_25), ...
        n_nb, n_g);
end

fprintf('\nDone.\n');

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function peak_pos = method_robust_peak(A)
% M2: Position of maximum amplitude
    [~, peak_pos] = max(A);
end

function [peak_pos, no_bump] = method_flanked_peak(A)
% M3: Flanked peak — peak must have top-4 ranked neighbors on both sides.
%
%   The peak (rank 1 position) must:
%     (a) not be on either edge (pos 1 or 11)
%     (b) have at least one of rank 2/3/4 positions as its LEFT neighbor
%     (c) have at least one of rank 2/3/4 positions as its RIGHT neighbor
%
%   If conditions fail, returns NaN and no_bump = true.

    n_pos = numel(A);
    [~, sorted_idx] = sort(A, 'descend');
    peak = sorted_idx(1);  % position with rank 1

    % Check: peak not on edge
    if peak == 1 || peak == n_pos
        peak_pos = NaN;
        no_bump = true;
        return;
    end

    % Get ranks of left and right neighbors
    ranks = zeros(1, n_pos);
    for r = 1:n_pos
        ranks(sorted_idx(r)) = r;
    end

    left_rank  = ranks(peak - 1);
    right_rank = ranks(peak + 1);

    % Both neighbors must be in top 4
    if left_rank <= 4 && right_rank <= 4
        peak_pos = peak;
        no_bump = false;
    else
        peak_pos = NaN;
        no_bump = true;
    end
end

function centroid = method_full_centroid(A)
% M4: Response-weighted centroid over all depolarizing positions (A > 0).
    pos = 1:numel(A);
    mask = A > 0;
    if any(mask)
        centroid = sum(A(mask) .* pos(mask)) / sum(A(mask));
    else
        centroid = 6;  % fallback
    end
end

function centroid = method_fwhm_centroid(A)
% M5: FWHM bump centroid — find contiguous half-max region around peak,
%     compute center-of-mass over those positions only.
    n_pos = numel(A);
    [peak_val, peak_pos] = max(A);

    if peak_val <= 0
        centroid = 6;  % fallback
        return;
    end

    threshold = peak_val / 2;

    % Grow left
    left = peak_pos;
    while left > 1 && A(left - 1) >= threshold
        left = left - 1;
    end

    % Grow right
    right = peak_pos;
    while right < n_pos && A(right + 1) >= threshold
        right = right + 1;
    end

    bump = left:right;
    A_bump = A(bump);
    centroid = sum(A_bump .* bump) / sum(A_bump);
end

function result = ternary(condition, val_true, val_false)
% Helper: ternary operator
    if condition
        result = val_true;
    else
        result = val_false;
    end
end
