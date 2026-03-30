% DIAGNOSE_OFF_DIRECTION_SWAP  Compare LUT direction mappings between ON
% and OFF cells to find the source of the swapped direction pairs.
%
%   For a sample ON cell and all OFF cells, prints:
%   - pattern_order and func_order (full and slow-bar-filtered)
%   - The slow bar mask (which indices are selected)
%   - Data row → LUT direction mapping
%   - The actual segment durations from the parser (to verify correct extraction)

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

%% Pick sample cells: 1 ON control, 1 OFF control, 1 OFF TTL
test_cells = {
    fullfile(data_root, 'control', 'ON',  '2025_07_25_16_34'), 'ON ctrl';
    fullfile(data_root, 'control', 'OFF', '2025_07_17_14_38'), 'OFF ctrl';
    fullfile(data_root, 'ttl',     'OFF', '2025_07_29_14_11'), 'OFF TTL';
};

orig_dir = pwd;

for c = 1:size(test_cells, 1)
    exp_folder = test_cells{c, 1};
    label = test_cells{c, 2};

    fprintf('\n%s\n', repmat('=', 1, 72));
    fprintf('CELL: %s  (%s)\n', label, exp_folder);
    fprintf('%s\n\n', repmat('=', 1, 72));

    % Load metadata
    ce = load(fullfile(exp_folder, 'currentExp.mat'), ...
        'pattern_order', 'func_order', 'metadata');

    fprintf('Strain: %s\n', ce.metadata.Strain);
    fprintf('Frame:  %d\n', ce.metadata.Frame);
    fprintf('pattern_order (%d entries): %s\n', numel(ce.pattern_order), mat2str(ce.pattern_order));
    fprintf('func_order    (%d entries): %s\n', numel(ce.func_order), mat2str(ce.func_order));

    % Show the slow bar mask
    slow_func_mask = (ce.func_order == 3 | ce.func_order == 4);
    slow_indices = find(slow_func_mask);
    fprintf('\nSlow bar mask (func==3|4) selects indices: %s\n', mat2str(slow_indices));
    fprintf('  -> %d slow bar conditions\n', sum(slow_func_mask));

    slow_patterns = ce.pattern_order(slow_func_mask);
    slow_funcs = ce.func_order(slow_func_mask);

    fprintf('\nSlow bar (pattern, func) pairs in metadata order:\n');
    for i = 1:numel(slow_patterns)
        row_mask = (Tbl.pattern == slow_patterns(i)) & (Tbl.function == slow_funcs(i));
        if any(row_mask)
            dir_deg = Tbl.direction(row_mask);
            orient_deg = Tbl.orientation(row_mask);
        else
            dir_deg = NaN;
            orient_deg = NaN;
        end
        fprintf('  DataRow %2d: pat=%2d func=%d -> dir=%5.1f° orient=%5.1f°\n', ...
            i, slow_patterns(i), slow_funcs(i), dir_deg, orient_deg);
    end

    % Now show what plot_order maps to
    fprintf('\nSubplot position -> DataRow -> LUT direction:\n');
    [lut_directions, ~, ~, ~] = ...
        verify_lut_directions(Tbl, ce.pattern_order, ce.func_order, plot_order);
    lut_dirs_ordered = lut_directions(plot_order);
    for sp = 1:16
        dr = plot_order(sp);
        fprintf('  Subplot %2d -> DataRow %2d -> %5.1f°\n', sp, dr, lut_dirs_ordered(sp));
    end

    % Check for direction ordering: are directions monotonically spaced?
    fprintf('\nDirections by data row (1-16): ');
    fprintf('%.1f ', lut_directions);
    fprintf('\n');

    % Load actual recording and check segment extraction
    fprintf('\nParsing actual recording to verify segment order...\n');
    cd(exp_folder);
    [~, ~, Log, ~, ~] = load_protocol2_data(exp_folder);
    cd(orig_dir);

    f_data = Log.ADC.Volts(1, :);
    v_data = Log.ADC.Volts(2, :) * 10;

    % Quick segment analysis for cycle 1 only
    % Find cycle boundaries
    zero_mask = (f_data == 0);
    d_zero = diff([0 zero_mask 0]);
    gap_starts = find(d_zero == 1);
    gap_ends = find(d_zero == -1) - 1;
    gap_lens = gap_ends - gap_starts + 1;
    long_gaps = find(gap_lens >= 80000);

    if numel(long_gaps) >= 2
        c1_start = gap_ends(long_gaps(1)) + 1;
        c1_end = gap_starts(long_gaps(2)) - 1;
        fprintf('Cycle 1: samples %d to %d (%.1f s)\n', c1_start, c1_end, (c1_end-c1_start+1)/10000);

        f_cycle = f_data(c1_start:c1_end);

        % Find transitions
        df = abs(diff(f_cycle));
        trans_local = find(df > 9);
        trans_abs = c1_start + trans_local;

        % Build segments
        seg_starts = [c1_start, trans_abs + 1];  % +1: first sample after transition
        seg_ends = [trans_abs, c1_end];

        n_segs = numel(seg_starts);
        fprintf('Found %d segments in cycle 1\n', n_segs);

        fprintf('\nSegment analysis (cycle 1):\n');
        fprintf('  %4s  %10s  %8s  %8s  %6s  %5s\n', ...
            'Seg', 'Duration', 'FrameStd', 'MeanFrm', 'Moving', 'Slow');
        slow_count = 0;
        for s = 1:n_segs
            seg_f = f_data(seg_starts(s):seg_ends(s));
            dur = numel(seg_f);
            f_std = std(double(seg_f));
            f_mean = mean(double(seg_f));
            is_moving = f_std > 5;
            is_slow = dur > 15000;
            if is_moving && is_slow
                slow_count = slow_count + 1;
                marker = sprintf('<-- SLOW #%d', slow_count);
            elseif is_moving
                marker = '    (fast)';
            else
                marker = '';
            end
            fprintf('  %4d  %10d  %8.1f  %8.1f  %6s  %5s  %s\n', ...
                s, dur, f_std, f_mean, string(is_moving), string(is_slow), marker);
        end
        fprintf('Total slow-moving segments in cycle 1: %d\n', slow_count);
    end

    fprintf('\n');
end

fprintf('Done.\n');
