% COMPARE_M6_SLOW_VS_FAST  Compare M6 RF center between slow and fast flash speeds.
%
%   For each of the 25 late-batch cells, extracts both slow (80ms) and fast
%   (14ms) bar flash traces, computes the M6 (68%-area) centroid independently
%   for each speed, and evaluates three alignment strategies:
%     1. Slow-M6 for both speeds
%     2. Fast-M6 for both speeds
%     3. Combined-M6 (average of normalized amplitude profiles)
%
%   Usage:
%     addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
%     run('scripts/compare_m6_slow_vs_fast.m');

%% Setup
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root = '/Users/reiserm/Documents/ttl_1DRF';

%% Load batch results (for pd_info per cell)
res_file = fullfile(data_root, 'population_results', 'batch_results.mat');
fprintf('Loading batch results: %s\n', res_file);
S = load(res_file, 'results');
results = S.results;
n_cells = numel(results);
fprintf('Loaded %d cells.\n\n', n_cells);

%% Timing constants
% Slow flash (80ms)
SLOW.bl_samples  = 1:5000;
SLOW.resp_start  = 5001;
SLOW.resp_end    = 6551;   % onset + 1550
SLOW.pctile      = 99.5;

% Fast flash (14ms) — prop_int=0.5 → gap=2500
FAST.bl_samples  = 1:2500;
FAST.resp_start  = 2501;
FAST.resp_end    = 4051;   % onset + 1550
FAST.pctile      = 99.5;

AREA_FRAC = 0.68;  % M6 area fraction
prop_int  = 0.5;

%% Load LUT (needed for pd_info reconstruction)
script_dir = fileparts(mfilename('fullpath'));
lut_path = fullfile(fileparts(script_dir), 'src', 'analysis', 'protocol2', 'bar_lut.mat');

%% Process each cell
T = struct();  % results table

for ci = 1:n_cells
    r = results(ci);
    exp_folder = fullfile(data_root, r.folder);

    fprintf('[%2d/%d] %s ... ', ci, n_cells, r.date_str);

    % --- Load raw data (same as batch pipeline) ---
    orig_dir = pwd;
    [~, ~, Log, ~, ~] = load_protocol2_data(exp_folder);
    cd(orig_dir);
    f_data = Log.ADC.Volts(1, :);
    v_data = Log.ADC.Volts(2, :) * 10;

    % --- Parse both flash speeds ---
    [~, ~, mean_slow, mean_fast] = parse_bar_flash_data(f_data, v_data, prop_int);

    % --- Reconstruct pd_info (need bar_flash_col, pos_order) ---
    LUT = load(lut_path);
    plot_order = [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16];
    pattern_offset = 2;

    ce = load(fullfile(exp_folder, 'currentExp.mat'), 'pattern_order', 'func_order');
    Tbl_local = LUT.Tbl;
    [lut_directions, lut_orientations, lut_patterns, lut_functions] = ...
        verify_lut_directions(Tbl_local, ce.pattern_order, ce.func_order, plot_order);

    bar_data = parse_bar_data(f_data, v_data);
    sweep_opts.baseline_range = [1000 9000];
    sweep_opts.stim_trim_end  = 7000;
    sweep_opts.percentile     = 98;
    max_v = compute_bar_sweep_responses(bar_data, plot_order, sweep_opts);
    pd_info = find_pd_from_lut(max_v, lut_directions, lut_orientations, ...
        lut_patterns, lut_functions, plot_order, Tbl_local, pattern_offset);

    % --- Extract amplitude profiles ---
    A_slow = extract_amplitude_profile(mean_slow, pd_info.bar_flash_col, ...
        pd_info.pos_order, SLOW);
    A_fast = extract_amplitude_profile(mean_fast, pd_info.bar_flash_col, ...
        pd_info.pos_order, FAST);

    % Also do orthogonal axis
    A_slow_ortho = extract_amplitude_profile(mean_slow, pd_info.ortho_flash_col, ...
        pd_info.ortho_pos_order, SLOW);
    A_fast_ortho = extract_amplitude_profile(mean_fast, pd_info.ortho_flash_col, ...
        pd_info.ortho_pos_order, FAST);

    % --- Compute M6 for each speed independently ---
    m6_slow = compute_m6_centroid(max(A_slow, 0), AREA_FRAC);
    m6_fast = compute_m6_centroid(max(A_fast, 0), AREA_FRAC);

    % --- Combined M6: normalize to peak=1, average, then M6 ---
    A_slow_norm = A_slow / max(max(A_slow, 0));
    A_fast_norm = A_fast / max(max(A_fast, 0));
    % Handle edge case: if peak is 0 or negative
    if max(A_slow) <= 0, A_slow_norm = zeros(size(A_slow)); end
    if max(A_fast) <= 0, A_fast_norm = zeros(size(A_fast)); end
    A_combined = (A_slow_norm + A_fast_norm) / 2;
    m6_comb = compute_m6_centroid(max(A_combined, 0), AREA_FRAC);

    % --- Ortho M6 ---
    m6_slow_ortho = compute_m6_centroid(max(A_slow_ortho, 0), AREA_FRAC);
    m6_fast_ortho = compute_m6_centroid(max(A_fast_ortho, 0), AREA_FRAC);

    % --- Store ---
    T(ci).date_str   = r.date_str;
    T(ci).is_on      = r.is_on;
    T(ci).is_ttl     = r.is_ttl;
    T(ci).group      = get_group_label(r.is_on, r.is_ttl);

    % PD axis
    T(ci).m6_slow_centroid   = m6_slow.centroid;
    T(ci).m6_slow_int        = m6_slow.centroid_int;
    T(ci).m6_slow_bw         = m6_slow.bump_width;
    T(ci).m6_fast_centroid   = m6_fast.centroid;
    T(ci).m6_fast_int        = m6_fast.centroid_int;
    T(ci).m6_fast_bw         = m6_fast.bump_width;
    T(ci).m6_comb_centroid   = m6_comb.centroid;
    T(ci).m6_comb_int        = m6_comb.centroid_int;
    T(ci).m6_comb_bw         = m6_comb.bump_width;

    T(ci).delta_slow_fast    = abs(m6_slow.centroid - m6_fast.centroid);
    T(ci).delta_int          = abs(m6_slow.centroid_int - m6_fast.centroid_int);

    % Ortho axis
    T(ci).m6_slow_ortho_int  = m6_slow_ortho.centroid_int;
    T(ci).m6_fast_ortho_int  = m6_fast_ortho.centroid_int;
    T(ci).delta_ortho_int    = abs(m6_slow_ortho.centroid_int - m6_fast_ortho.centroid_int);

    % Amplitude comparison
    T(ci).peak_slow = max(A_slow);
    T(ci).peak_fast = max(A_fast);
    T(ci).ratio_fast_slow = max(A_fast) / max(max(A_slow, 0.01));

    % Raw profiles for later inspection
    T(ci).A_slow = A_slow;
    T(ci).A_fast = A_fast;

    fprintf('slow_M6=%.1f(%d) fast_M6=%.1f(%d) comb_M6=%.1f(%d) delta=%.1f\n', ...
        m6_slow.centroid, m6_slow.centroid_int, ...
        m6_fast.centroid, m6_fast.centroid_int, ...
        m6_comb.centroid, m6_comb.centroid_int, ...
        T(ci).delta_slow_fast);
end

%% ========================= Summary Table ================================
fprintf('\n\n');
fprintf('=================================================================================================================================\n');
fprintf('%-25s %-10s  Slow_M6  Fast_M6  Comb_M6  SlowInt FastInt CombInt  Delta  DeltaInt  SlowBW FastBW  PkSlow  PkFast  Ratio\n', 'Cell', 'Group');
fprintf('=================================================================================================================================\n');

for ci = 1:n_cells
    fprintf('%-25s %-10s  %5.2f    %5.2f    %5.2f      %2d      %2d      %2d    %5.2f    %2d      %2d     %2d    %6.1f   %6.1f   %.2f\n', ...
        T(ci).date_str, T(ci).group, ...
        T(ci).m6_slow_centroid, T(ci).m6_fast_centroid, T(ci).m6_comb_centroid, ...
        T(ci).m6_slow_int, T(ci).m6_fast_int, T(ci).m6_comb_int, ...
        T(ci).delta_slow_fast, T(ci).delta_int, ...
        T(ci).m6_slow_bw, T(ci).m6_fast_bw, ...
        T(ci).peak_slow, T(ci).peak_fast, T(ci).ratio_fast_slow);
end

%% ========================= Summary Statistics ===========================
fprintf('\n\n=== PD Axis: M6 Agreement Summary ===\n');
deltas     = [T.delta_slow_fast];
delta_ints = [T.delta_int];
fprintf('  Fractional centroid delta (slow vs fast):\n');
fprintf('    Mean: %.2f positions\n', mean(deltas));
fprintf('    Median: %.2f positions\n', median(deltas));
fprintf('    Max: %.2f positions\n', max(deltas));
fprintf('  Integer centroid agreement:\n');
fprintf('    Exact match: %d / %d (%.0f%%)\n', sum(delta_ints == 0), n_cells, 100*mean(delta_ints == 0));
fprintf('    Off by 1: %d / %d\n', sum(delta_ints == 1), n_cells);
fprintf('    Off by >=2: %d / %d\n', sum(delta_ints >= 2), n_cells);

fprintf('\n=== Ortho Axis: M6 Agreement Summary ===\n');
delta_ortho = [T.delta_ortho_int];
fprintf('  Integer centroid agreement:\n');
fprintf('    Exact match: %d / %d (%.0f%%)\n', sum(delta_ortho == 0), n_cells, 100*mean(delta_ortho == 0));
fprintf('    Off by 1: %d / %d\n', sum(delta_ortho == 1), n_cells);
fprintf('    Off by >=2: %d / %d\n', sum(delta_ortho >= 2), n_cells);

fprintf('\n=== Amplitude Comparison ===\n');
fprintf('  Mean peak slow: %.1f mV\n', mean([T.peak_slow]));
fprintf('  Mean peak fast: %.1f mV\n', mean([T.peak_fast]));
fprintf('  Mean ratio (fast/slow): %.2f\n', mean([T.ratio_fast_slow]));

fprintf('\n=== Bump Width Comparison ===\n');
fprintf('  Mean slow BW: %.1f positions\n', mean([T.m6_slow_bw]));
fprintf('  Mean fast BW: %.1f positions\n', mean([T.m6_fast_bw]));

%% ========================= By Group =====================================
groups = {'ON ctrl', 'ON TTL', 'OFF ctrl', 'OFF TTL'};
masks = {[T.is_on] & ~[T.is_ttl], [T.is_on] & [T.is_ttl], ...
         ~[T.is_on] & ~[T.is_ttl], ~[T.is_on] & [T.is_ttl]};

fprintf('\n=== By Group ===\n');
fprintf('%-10s  n   Delta(mean)  Delta(med)  Match%%  SlowBW  FastBW  PkRatio\n', 'Group');
for g = 1:4
    m = masks{g};
    if sum(m) == 0, continue; end
    d = deltas(m);
    di = delta_ints(m);
    fprintf('%-10s %2d   %5.2f        %5.2f       %3.0f%%    %4.1f    %4.1f    %.2f\n', ...
        groups{g}, sum(m), mean(d), median(d), 100*mean(di==0), ...
        mean([T(m).m6_slow_bw]), mean([T(m).m6_fast_bw]), mean([T(m).ratio_fast_slow]));
end

%% ========================= Alignment Strategy Comparison ================
% For each strategy, compute the "alignment error" = how far the chosen
% center is from both speeds' true M6. Lower is better.
fprintf('\n\n=== Alignment Strategy Evaluation ===\n');
fprintf('For each strategy, report mean distance from chosen center to each speed''s M6:\n\n');

slow_c = [T.m6_slow_centroid];
fast_c = [T.m6_fast_centroid];
comb_c = [T.m6_comb_centroid];

% Strategy 1: use slow M6
err1_to_slow = mean(abs(slow_c - slow_c));  % = 0 by definition
err1_to_fast = mean(abs(slow_c - fast_c));
fprintf('  Strategy 1 (Slow-M6 for both):\n');
fprintf('    Mean dist to slow M6: %.3f (exact by definition)\n', err1_to_slow);
fprintf('    Mean dist to fast M6: %.3f\n', err1_to_fast);

% Strategy 2: use fast M6
err2_to_slow = mean(abs(fast_c - slow_c));
err2_to_fast = mean(abs(fast_c - fast_c));  % = 0
fprintf('  Strategy 2 (Fast-M6 for both):\n');
fprintf('    Mean dist to slow M6: %.3f\n', err2_to_slow);
fprintf('    Mean dist to fast M6: %.3f (exact by definition)\n', err2_to_fast);

% Strategy 3: use combined M6
err3_to_slow = mean(abs(comb_c - slow_c));
err3_to_fast = mean(abs(comb_c - fast_c));
fprintf('  Strategy 3 (Combined-M6 for both):\n');
fprintf('    Mean dist to slow M6: %.3f\n', err3_to_slow);
fprintf('    Mean dist to fast M6: %.3f\n', err3_to_fast);
fprintf('    Mean total dist: %.3f\n', (err3_to_slow + err3_to_fast)/2);

% Recommendation
fprintf('\n=== Recommendation ===\n');
if max(deltas) < 1.0 && mean(delta_ints == 0) > 0.7
    fprintf('  M6 is highly consistent between speeds (max delta %.2f, %.0f%% exact match).\n', ...
        max(deltas), 100*mean(delta_ints == 0));
    fprintf('  --> Use SLOW-M6 for both speeds (Strategy 1).\n');
    fprintf('      Rationale: slow flash has higher SNR, and M6 centers agree well.\n');
elseif mean(deltas) < 1.5
    fprintf('  M6 is moderately consistent (mean delta %.2f).\n', mean(deltas));
    fprintf('  --> Use COMBINED-M6 (Strategy 3) as a compromise.\n');
else
    fprintf('  M6 shows substantial disagreement (mean delta %.2f).\n', mean(deltas));
    fprintf('  --> Use COMBINED-M6 (Strategy 3) to split the difference,\n');
    fprintf('      or consider independent alignment per speed.\n');
end

fprintf('\nDone.\n');


%% ========================= Local Functions ==============================

function A = extract_amplitude_profile(mean_bf, flash_col, pos_order, timing)
% Extract 1x11 amplitude profile from mean bar flash data.
    n_pos = 11;
    A = zeros(1, n_pos);
    for pos_idx = 1:n_pos
        flash_pos = pos_order(pos_idx);
        ts = mean_bf{flash_pos, flash_col};
        if ~isempty(ts)
            bl_end = min(max(timing.bl_samples), numel(ts));
            bl_mean = mean(ts(1:bl_end));
            win_end = min(timing.resp_end, numel(ts));
            resp_win = timing.resp_start : win_end;
            if ~isempty(resp_win)
                A(pos_idx) = prctile(ts(resp_win) - bl_mean, timing.pctile);
            end
        end
    end
end


function m6 = compute_m6_centroid(A, area_fraction)
% M6: 68%-area bump centroid via greedy expansion from peak.
    n_pos = numel(A);
    [peak_val, peak_pos] = max(A);
    total_area = sum(A);

    if peak_val <= 0 || total_area <= 0
        m6.centroid     = 6;
        m6.centroid_int = 6;
        m6.bump_range   = [1, n_pos];
        m6.bump_width   = n_pos;
        return;
    end

    target = area_fraction * total_area;
    left  = peak_pos;
    right = peak_pos;
    cum   = A(peak_pos);

    while cum < target && (left > 1 || right < n_pos)
        can_left  = (left > 1);
        can_right = (right < n_pos);
        if can_left && can_right
            if A(left - 1) >= A(right + 1)
                left = left - 1;
                cum  = cum + A(left);
            else
                right = right + 1;
                cum   = cum + A(right);
            end
        elseif can_left
            left = left - 1;
            cum  = cum + A(left);
        else
            right = right + 1;
            cum   = cum + A(right);
        end
    end

    bump   = left:right;
    A_bump = A(bump);
    m6.centroid     = sum(A_bump .* bump) / sum(A_bump);
    m6.centroid_int = round(m6.centroid);
    m6.centroid_int = max(1, min(n_pos, m6.centroid_int));
    m6.bump_range   = [left, right];
    m6.bump_width   = right - left + 1;
end


function lbl = get_group_label(is_on, is_ttl)
    if is_on && ~is_ttl
        lbl = 'ON ctrl';
    elseif is_on && is_ttl
        lbl = 'ON TTL';
    elseif ~is_on && ~is_ttl
        lbl = 'OFF ctrl';
    else
        lbl = 'OFF TTL';
    end
end
