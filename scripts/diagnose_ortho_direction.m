% DIAGNOSE_ORTHO_DIRECTION  Check orthogonal bar position ordering for all 25 cells.
%
%   For each cell, computes:
%     - PD direction (from vector sum)
%     - Orthogonal bar orientation and forward sweep direction (from LUT)
%     - Right-hand rule target direction (PD - 90°, i.e. 90° clockwise from PD)
%     - Whether current pos_order matches the right-hand rule
%
%   Also generates a diagnostic polar plot showing PD, ortho axis, and
%   right-hand direction for a few example cells.
%
%   USAGE:
%     run('scripts/diagnose_ortho_direction.m')
%
%   REQUIRES: src/ on path, bar_lut.mat, CircStat toolbox

clear; close all;
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root = '/Users/reiserm/Documents/ttl_1DRF';
lut_file  = fullfile(fileparts(mfilename('fullpath')), ...
    '..', 'src', 'analysis', 'protocol2', 'bar_lut.mat');
lut = load(lut_file, 'Tbl');
Tbl = lut.Tbl;

% Default options (from batch_analyze_1DRF)
plot_order = [13 14 15 16 9 10 11 12 5 6 7 8 1 2 3 4];
pattern_offset = 2;
baseline_range = [1000, 9000];
stim_trim_end  = 1000;
pctile_val     = 99.5;

% Get all experiment folders
exp_dirs = dir(fullfile(data_root, '2025_*'));
exp_dirs = exp_dirs([exp_dirs.isdir]);
fprintf('Found %d experiment folders\n\n', numel(exp_dirs));

% Print header
fprintf('%-20s  %6s  %6s  %8s  %8s  %8s  %8s  %5s  %s\n', ...
    'Experiment', 'PD_dir', 'PD_fn', 'Ortho_or', 'Ortho_fw', 'RH_targ', 'AngDist', 'Match', 'pos_order_needed');
fprintf('%s\n', repmat('-', 1, 110));

results_table = [];

for i = 1:numel(exp_dirs)
    exp_folder = fullfile(data_root, exp_dirs(i).name);

    try
        % Load data
        [date_str, ~, Log, ~, ~] = load_protocol2_data(exp_folder);
        f_data = Log.ADC.Volts(1, :);
        v_data = Log.ADC.Volts(2, :) * 10;

        ce = load(fullfile(exp_folder, 'currentExp.mat'), ...
            'pattern_order', 'func_order', 'metadata');

        % Verify LUT
        [lut_directions, lut_orientations, lut_patterns, lut_functions] = ...
            verify_lut_directions(Tbl, ce.pattern_order, ce.func_order, plot_order);

        % Compute bar sweep responses
        bar_data = parse_bar_data(f_data, v_data);
        sweep_opts.baseline_range = baseline_range;
        sweep_opts.stim_trim_end  = stim_trim_end;
        sweep_opts.percentile     = pctile_val;
        max_v = compute_bar_sweep_responses(bar_data, plot_order, sweep_opts);

        % Find PD (this gives us pd_info with all the LUT lookups)
        pd_info = find_pd_from_lut(max_v, lut_directions, lut_orientations, ...
            lut_patterns, lut_functions, plot_order, Tbl, pattern_offset);

        % --- Ortho direction analysis ---
        % Look up orthogonal bar's forward sweep direction from LUT
        ortho_exp_pat = pd_info.ortho_flash_col + pattern_offset;
        ortho_fwd_mask = Tbl.pattern == ortho_exp_pat & Tbl.function == 3;  % function 3 = forward
        ortho_rev_mask = Tbl.pattern == ortho_exp_pat & Tbl.function == 4;  % function 4 = reverse

        if any(ortho_fwd_mask)
            ortho_dir_fwd = Tbl.direction(ortho_fwd_mask);
        else
            % Try odd functions more broadly
            ortho_odd_mask = Tbl.pattern == ortho_exp_pat & mod(Tbl.function, 2) == 1;
            if any(ortho_odd_mask)
                ortho_dir_fwd = Tbl.direction(find(ortho_odd_mask, 1));
            else
                fprintf('%-20s  ** No forward function found for ortho pattern %d **\n', ...
                    exp_dirs(i).name, ortho_exp_pat);
                continue;
            end
        end

        % Right-hand rule: 90° clockwise from PD
        rh_target = mod(pd_info.pd_direction - 90, 360);

        % Angular distance between ortho forward direction and RH target
        ang_diff = abs(ortho_dir_fwd - rh_target);
        ang_diff = min(ang_diff, 360 - ang_diff);  % handle wraparound

        % If ang_diff < 90°, forward matches RH rule → pos_order = 1:11
        % If ang_diff > 90°, reverse matches RH rule → pos_order = 11:-1:1
        if ang_diff < 90
            ortho_pos_order_needed = '1:11  (fwd)';
            rh_matches_fwd = true;
        else
            ortho_pos_order_needed = '11:-1:1 (rev)';
            rh_matches_fwd = false;
        end

        % Current pos_order (from PD function)
        is_pd_forward = mod(pd_info.pd_function, 2) == 1;
        if is_pd_forward
            current_matches = rh_matches_fwd;  % both forward
        else
            current_matches = ~rh_matches_fwd;  % both reverse
        end

        fprintf('%-20s  %6.1f  %6d  %8.1f  %8.1f  %8.1f  %8.1f  %5s  %s\n', ...
            exp_dirs(i).name, pd_info.pd_direction, pd_info.pd_function, ...
            pd_info.ortho_orientation, ortho_dir_fwd, rh_target, ang_diff, ...
            yesno(current_matches), ortho_pos_order_needed);

        % Store for summary
        row.name = exp_dirs(i).name;
        row.pd_direction = pd_info.pd_direction;
        row.pd_function = pd_info.pd_function;
        row.ortho_orientation = pd_info.ortho_orientation;
        row.ortho_dir_fwd = ortho_dir_fwd;
        row.rh_target = rh_target;
        row.ang_diff = ang_diff;
        row.rh_matches_fwd = rh_matches_fwd;
        row.current_matches = current_matches;
        row.is_pd_forward = is_pd_forward;

        if isempty(results_table)
            results_table = row;
        else
            results_table(end+1) = row;  %#ok<SAGROW>
        end

    catch ME
        fprintf('%-20s  ** ERROR: %s **\n', exp_dirs(i).name, ME.message);
    end
end

% --- Summary ---
fprintf('\n=== SUMMARY ===\n');
n_total = numel(results_table);
n_match = sum([results_table.current_matches]);
n_mismatch = n_total - n_match;
fprintf('Total cells: %d\n', n_total);
fprintf('Current pos_order matches RH rule: %d (%.0f%%)\n', n_match, 100*n_match/n_total);
fprintf('Current pos_order MISMATCHES RH rule: %d (%.0f%%)\n', n_mismatch, 100*n_mismatch/n_total);

if n_mismatch > 0
    fprintf('\nMismatched cells (ortho positions currently flipped):\n');
    for k = 1:n_total
        if ~results_table(k).current_matches
            fprintf('  %s — PD: %.1f°, fn: %d, ortho_fwd: %.1f°, RH_target: %.1f°\n', ...
                results_table(k).name, results_table(k).pd_direction, ...
                results_table(k).pd_function, results_table(k).ortho_dir_fwd, ...
                results_table(k).rh_target);
        end
    end
end

% --- Diagnostic polar plot for 4 example cells ---
fig = figure('Name', 'Ortho Direction Diagnostic', 'Position', [50 50 1200 300]);
n_examples = min(4, n_total);
tiledlayout(1, n_examples, 'TileSpacing', 'compact', 'Padding', 'compact');

for k = 1:n_examples
    nexttile;
    r = results_table(k);

    % Polar plot showing PD and ortho directions
    pd_rad = deg2rad(r.pd_direction);
    rh_rad = deg2rad(r.rh_target);
    ortho_fwd_rad = deg2rad(r.ortho_dir_fwd);
    ortho_rev_rad = ortho_fwd_rad + pi;

    polarplot([pd_rad pd_rad], [0 1], 'r-', 'LineWidth', 2.5); hold on;
    polarplot([rh_rad rh_rad], [0 0.8], 'b-', 'LineWidth', 2);
    polarplot([ortho_fwd_rad ortho_rev_rad], [0.6 0.6], ...
        'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);

    % Arrow markers
    polarplot(pd_rad, 1, 'r^', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    polarplot(rh_rad, 0.8, 'bv', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

    % Mark which end is "pos 1→11" direction for ortho
    if r.rh_matches_fwd
        polarplot(ortho_fwd_rad, 0.6, 'g>', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
    else
        polarplot(ortho_rev_rad, 0.6, 'g>', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
    end

    title(sprintf('%s\nPD=%.0f° fn=%d %s', ...
        strrep(r.name, '_', '-'), r.pd_direction, r.pd_function, ...
        ternary(r.current_matches, 'OK', 'FLIP')), 'FontSize', 9);
    rlim([0 1.2]);
end

% Save diagnostic figure
preview_dir = fullfile(data_root, 'figure_previews');
if ~exist(preview_dir, 'dir'), mkdir(preview_dir); end
exportgraphics(fig, fullfile(preview_dir, 'diag_ortho_direction.png'), 'Resolution', 150);
fprintf('\nDiagnostic figure saved: %s\n', fullfile(preview_dir, 'diag_ortho_direction.png'));


function s = yesno(val)
    if val, s = 'YES'; else, s = 'NO'; end
end

function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end
