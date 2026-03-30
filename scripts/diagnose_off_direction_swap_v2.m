% DIAGNOSE_OFF_DIRECTION_SWAP_V2  Decode the actual bar direction for each
% slow segment by analyzing the frame signal sweep direction.
%
%   For each slow segment, computes:
%     - Frame value at start vs end (indicates sweep direction)
%     - Mean frame during first/last 500 samples
%   Compares ON vs OFF cells to find where directions diverge.

%% Setup
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root = '/Users/reiserm/Documents/ttl_1DRF/pre-bar-flash';
plot_order = [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16];

% Load LUT
lut_path = fullfile(fileparts(mfilename('fullpath')), ...
    '..', 'src', 'analysis', 'protocol2', 'bar_lut.mat');
lut_path = char(java.io.File(lut_path).getCanonicalPath());
S_lut = load(lut_path, 'Tbl');
Tbl = S_lut.Tbl;

%% Test cells
test_cells = {
    fullfile(data_root, 'control', 'ON',  '2025_07_25_16_34'), 'ON ctrl';
    fullfile(data_root, 'control', 'ON',  '2025_07_22_16_38'), 'ON ctrl 2';
    fullfile(data_root, 'control', 'OFF', '2025_07_17_14_38'), 'OFF ctrl';
    fullfile(data_root, 'control', 'OFF', '2025_07_22_18_15'), 'OFF ctrl 2';
    fullfile(data_root, 'ttl',     'OFF', '2025_07_29_14_11'), 'OFF TTL';
};

orig_dir = pwd;

for c = 1:size(test_cells, 1)
    exp_folder = test_cells{c, 1};
    label = test_cells{c, 2};

    fprintf('\n%s\n', repmat('=', 1, 80));
    fprintf('CELL: %s\n', label);
    fprintf('%s\n', repmat('=', 1, 80));

    % Load data
    cd(exp_folder);
    [~, ~, Log, ~, ~] = load_protocol2_data(exp_folder);
    cd(orig_dir);

    f_data = Log.ADC.Volts(1, :);
    v_data = Log.ADC.Volts(2, :) * 10;

    ce = load(fullfile(exp_folder, 'currentExp.mat'), ...
        'pattern_order', 'func_order', 'metadata');
    fprintf('Frame: %d\n', ce.metadata.Frame);

    % Get LUT mapping
    [lut_directions, ~, lut_patterns, lut_functions] = ...
        verify_lut_directions(Tbl, ce.pattern_order, ce.func_order, plot_order);

    % Parse: find cycle boundaries and extract slow segments from cycle 1
    zero_mask = (f_data == 0);
    d_zero = diff([0 zero_mask 0]);
    gap_starts = find(d_zero == 1);
    gap_ends = find(d_zero == -1) - 1;
    gap_lens = gap_ends - gap_starts + 1;
    long_gaps = find(gap_lens >= 80000);

    c1_start = gap_ends(long_gaps(1)) + 1;
    c1_end = gap_starts(long_gaps(2)) - 1;

    f_cycle = f_data(c1_start:c1_end);

    % Find transitions and segments
    df = abs(diff(f_cycle));
    trans_local = find(df > 9);
    trans_abs = c1_start + trans_local;

    seg_starts = [c1_start, trans_abs + 1];
    seg_ends = [trans_abs, c1_end];

    % Classify and extract slow moving segments
    slow_segs = [];
    for s = 1:numel(seg_starts)
        seg_f = f_data(seg_starts(s):seg_ends(s));
        dur = numel(seg_f);
        f_std = std(double(seg_f));
        is_moving = f_std > 5;
        is_slow = dur > 15000;
        if is_moving && is_slow
            slow_segs(end+1,:) = [seg_starts(s), seg_ends(s)]; %#ok<SAGROW>
        end
    end

    fprintf('Found %d slow segments in cycle 1\n\n', size(slow_segs, 1));

    % For each slow segment, analyze the frame signal
    fprintf('%4s  %8s  %8s  %8s  %6s  %8s  %8s  %8s\n', ...
        'Seg#', 'F_start', 'F_end', 'F_delta', 'F_dir', ...
        'LUT_dir', 'LUT_pat', 'LUT_fun');
    fprintf('%s\n', repmat('-', 1, 72));

    n_edge = 500;  % samples to average at start/end

    for i = 1:min(16, size(slow_segs, 1))
        s_start = slow_segs(i, 1);
        s_end = slow_segs(i, 2);
        seg_f = f_data(s_start:s_end);

        f_start = mean(seg_f(1:n_edge));
        f_end = mean(seg_f(end-n_edge+1:end));
        f_delta = f_end - f_start;

        if f_delta > 0
            f_dir_str = '+';
        else
            f_dir_str = '-';
        end

        fprintf('%4d  %8.2f  %8.2f  %8.2f  %6s  %8.1f  %8d  %8d\n', ...
            i, f_start, f_end, f_delta, f_dir_str, ...
            lut_directions(i), lut_patterns(i), lut_functions(i));
    end

    % Now print it organized by pairs (odd/even data rows)
    fprintf('\n--- Organized by pattern pairs ---\n');
    fprintf('%8s  %8s  %8s  |  %8s  %8s  %8s\n', ...
        'Pair', 'Seg(f3)', 'F_delta', 'Seg(f4)', 'F_delta', 'Pattern');
    fprintf('%s\n', repmat('-', 1, 60));

    for pair = 1:8
        idx_f3 = 2*pair - 1;  % odd rows = func 3
        idx_f4 = 2*pair;      % even rows = func 4

        sf3 = f_data(slow_segs(idx_f3,1):slow_segs(idx_f3,2));
        sf4 = f_data(slow_segs(idx_f4,1):slow_segs(idx_f4,2));

        delta_f3 = mean(sf3(end-n_edge+1:end)) - mean(sf3(1:n_edge));
        delta_f4 = mean(sf4(end-n_edge+1:end)) - mean(sf4(1:n_edge));

        fprintf('%8d  %8d  %+8.2f  |  %8d  %+8.2f  pat=%d\n', ...
            pair, idx_f3, delta_f3, idx_f4, delta_f4, lut_patterns(idx_f3));
    end

    fprintf('\n');
end

fprintf('Done.\n');
