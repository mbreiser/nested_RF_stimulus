% DIAGNOSE_PRE_BAR_FLASH  Validate stimulus structure for pre-bar-flash dataset.
%
%   Loads sample experiments from the pre-bar-flash dataset, detects grey-
%   screen gaps in the frame signal, prints gap structure, validates LUT
%   directions, and checks the pattern 0015/0016 issue for OFF cells.
%
%   Key finding from initial run: pre-bar-flash experiments have NO 3s
%   intra-cycle grey gaps. Each cycle runs continuously:
%     [10s grey] -> flashes -> 28dps bars -> 56dps bars -> [10s grey]
%   Only the 10s grey screens between cycles are frame==0 gaps.
%
%   See also PARSE_BAR_DATA, VERIFY_LUT_DIRECTIONS, BATCH_ANALYZE_PRE_BAR_FLASH

addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root = '/Users/reiserm/Documents/ttl_1DRF/pre-bar-flash';

% Load LUT
lut_path = fullfile(fileparts(which('batch_analyze_1DRF')), 'bar_lut.mat');
S_lut = load(lut_path, 'Tbl');
Tbl = S_lut.Tbl;

% Default plot_order (same as main pipeline)
plot_order = [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16];

% Sample experiments
sample_folders = {
    fullfile(data_root, 'control', 'ON',  '2025_07_22_16_38'), ...  % pre-Oct-15, ON ctrl
    fullfile(data_root, 'control', 'OFF', '2025_07_17_14_38'), ...  % pre-Oct-15, OFF ctrl
    fullfile(data_root, 'control', 'OFF', '2025_10_22_13_10')  ...  % post-Oct-15, OFF ctrl
};
sample_labels = {'Pre-Oct15 ON ctrl', 'Pre-Oct15 OFF ctrl', 'Post-Oct15 OFF ctrl'};

fprintf('========================================================================\n');
fprintf('  DIAGNOSTIC: Pre-Bar-Flash Dataset Stimulus Structure\n');
fprintf('  Date: %s\n', datestr(now));
fprintf('========================================================================\n');

%% 1. Count all experiments
fprintf('\n--- Experiment Count ---\n');
group_paths = {'control/ON', 'control/OFF', 'ttl/ON', 'ttl/OFF'};
for g = 1:numel(group_paths)
    d = dir(fullfile(data_root, group_paths{g}));
    d = d([d.isdir] & ~startsWith({d.name}, '.'));
    fprintf('  %s: %d experiments\n', group_paths{g}, numel(d));
end

%% 2. Analyze each sample experiment
for s = 1:numel(sample_folders)
    exp_folder = sample_folders{s};
    fprintf('\n========================================================================\n');
    fprintf('  SAMPLE %d: %s\n', s, sample_labels{s});
    fprintf('  Folder: %s\n', exp_folder);
    fprintf('========================================================================\n');

    % Load data
    orig_dir = pwd;
    [date_str, ~, Log, ~, ~] = load_protocol2_data(exp_folder);
    cd(orig_dir);

    f_data = Log.ADC.Volts(1, :);
    v_data = Log.ADC.Volts(2, :) * 10;

    % Load currentExp.mat
    ce = load(fullfile(exp_folder, 'currentExp.mat'));
    fprintf('\n  metadata.Strain: %s\n', ce.metadata.Strain);
    fprintf('  metadata.Frame:  %d\n', ce.metadata.Frame);
    fprintf('  Recording length: %.1f s (%d samples)\n', numel(f_data)/10000, numel(f_data));

    % Print pattern_order and func_order
    fprintf('\n  pattern_order (%d entries): [', numel(ce.pattern_order));
    fprintf('%d ', ce.pattern_order);
    fprintf(']\n');
    fprintf('  func_order    (%d entries): [', numel(ce.func_order));
    fprintf('%d ', ce.func_order);
    fprintf(']\n');

    % Count function files
    func_dir = fullfile(exp_folder, 'Functions');
    func_files = dir(fullfile(func_dir, '0*.mat'));
    fprintf('  Function files: %d\n', numel(func_files));
    for ff = 1:numel(func_files)
        fprintf('    %s\n', func_files(ff).name);
    end

    %% 2a. Detect ALL grey-screen gaps (frame == 0)
    fprintf('\n  --- Gap Detection (frame==0) ---\n');
    zero_mask = f_data == 0;
    d_mask = diff([0 zero_mask 0]);
    gap_starts = find(d_mask == 1);
    gap_ends = find(d_mask == -1) - 1;
    gap_lens = gap_ends - gap_starts + 1;

    fprintf('  All gaps where frame==0:\n');
    fprintf('  %-5s  %12s  %12s  %10s\n', 'Gap#', 'Start', 'End', 'Duration_s');
    fprintf('  %s\n', repmat('-', 1, 45));
    for g = 1:numel(gap_starts)
        dur_s = gap_lens(g) / 10000;
        if dur_s >= 0.5  % only show gaps > 0.5s
            fprintf('  %-5d  %12d  %12d  %10.2f\n', g, gap_starts(g), gap_ends(g), dur_s);
        end
    end

    % Long gaps (>= 3s)
    long_mask = gap_lens >= 30000;
    n_long_gaps = sum(long_mask);
    fprintf('\n  Gaps >= 3s: %d (these are cycle-delimiting grey screens)\n', n_long_gaps);

    %% 2b. Define cycle boundaries from long grey gaps
    long_gap_starts = gap_starts(long_mask);
    long_gap_ends = gap_ends(long_mask);
    long_gap_lens = gap_lens(long_mask);

    % Cycles run between the long grey screens
    % Cycle 1: from end of gap1 to start of gap2
    % Cycle 2: from end of gap2 to start of gap3
    % Cycle 3: from end of gap3 to start of gap4 (or end of recording)
    n_cycles = n_long_gaps - 1;  % last gap is trailing grey
    fprintf('  Number of stimulus cycles: %d\n', n_cycles);

    for c = 1:n_cycles
        cycle_start = long_gap_ends(c) + 1;
        cycle_end = long_gap_starts(c + 1) - 1;
        cycle_dur = (cycle_end - cycle_start + 1) / 10000;
        fprintf('  Cycle %d: samples %d to %d (%.1f s)\n', ...
            c, cycle_start, cycle_end, cycle_dur);
    end

    %% 2c. Analyze cycle 1 in detail — find transitions
    cycle1_start = long_gap_ends(1) + 1;
    cycle1_end = long_gap_starts(2) - 1;
    f_cycle = f_data(cycle1_start:cycle1_end);

    % Find ALL large transitions (abs(diff) > 9)
    transitions = find(abs(diff(f_cycle)) > 9);
    fprintf('\n  Cycle 1 analysis (total %.1f s):\n', numel(f_cycle)/10000);
    fprintf('  Total transitions (|diff|>9): %d\n', numel(transitions));

    % Time from cycle start to first transition
    if ~isempty(transitions)
        fprintf('  Time to first bar transition: %.2f s\n', transitions(1)/10000);
    end

    % Characterize segments between transitions
    % Build segments list
    segs = zeros(numel(transitions) + 1, 3); % [start, end, duration_ms]
    seg_starts = [1, transitions + 1];
    seg_ends = [transitions, numel(f_cycle)];
    n_segs = numel(seg_starts);

    % Classify segments by frame characteristics
    fprintf('\n  Segment analysis (first 70 of %d segments):\n', n_segs);
    fprintf('  %-5s  %10s  %10s  %10s  %12s  %12s  %10s\n', ...
        'Seg#', 'Start_s', 'End_s', 'Dur_ms', 'Mean_frame', 'Frame_std', 'Type');
    fprintf('  %s\n', repmat('-', 1, 75));

    bar_seg_count = 0;
    for seg_idx = 1:min(70, n_segs)
        seg_s = seg_starts(seg_idx);
        seg_e = seg_ends(seg_idx);
        seg_dur_ms = (seg_e - seg_s + 1) / 10;  % ms
        f_seg = f_cycle(seg_s:seg_e);
        f_mean = mean(f_seg);
        f_std = std(f_seg);

        % Classify:
        % - Flash: short duration (~160ms), small frame values
        % - Static bar: ~3s, steady frame value, std ~0
        % - Moving bar: ~1-2.3s, changing frame value, high std
        if f_std > 5
            seg_type = 'MOVING';
            bar_seg_count = bar_seg_count + 1;
        elseif seg_dur_ms > 2000
            seg_type = 'STATIC';
        elseif seg_dur_ms > 100
            seg_type = 'flash?';
        else
            seg_type = 'brief';
        end

        fprintf('  %-5d  %10.2f  %10.2f  %10.1f  %12.1f  %12.1f  %10s\n', ...
            seg_idx, seg_s/10000, seg_e/10000, seg_dur_ms, f_mean, f_std, seg_type);
    end

    fprintf('\n  Total MOVING segments in first 70: %d\n', bar_seg_count);

    %% 2d. Count all MOVING and STATIC bar segments across full cycle
    fprintf('\n  Full cycle segment classification:\n');
    moving_count = 0;
    static_count = 0;
    moving_durs = [];
    static_durs = [];

    for seg_idx = 1:n_segs
        seg_s = seg_starts(seg_idx);
        seg_e = seg_ends(seg_idx);
        seg_dur_ms = (seg_e - seg_s + 1) / 10;
        f_seg = f_cycle(seg_s:seg_e);
        f_std = std(f_seg);

        if f_std > 5 && seg_dur_ms > 500
            moving_count = moving_count + 1;
            moving_durs(end+1) = seg_dur_ms; %#ok<SAGROW>
        elseif f_std < 2 && seg_dur_ms > 2000
            static_count = static_count + 1;
            static_durs(end+1) = seg_dur_ms; %#ok<SAGROW>
        end
    end

    fprintf('  MOVING bar segments (std>5, dur>500ms): %d\n', moving_count);
    if ~isempty(moving_durs)
        fprintf('    Duration range: %.0f - %.0f ms\n', min(moving_durs), max(moving_durs));
        fprintf('    Mean duration: %.0f ms\n', mean(moving_durs));
    end
    fprintf('  STATIC bar segments (std<2, dur>2s): %d\n', static_count);
    if ~isempty(static_durs)
        fprintf('    Duration range: %.0f - %.0f ms\n', min(static_durs), max(static_durs));
    end

    % We expect: 16 slow bars (28dps) + 16 fast bars (56dps) = 32 moving segments
    % Each has a static interval before/after
    fprintf('  Expected: 32 moving segments (16 slow + 16 fast)\n');

    %% 2e. Try the parse_bar_data approach adapted for this protocol
    % Instead of using grey gaps as delimiters, find transitions in the FULL
    % recording and identify bar sweep regions by segment characteristics.
    %
    % Key insight: within each cycle, the 28dps bars come first (longer
    % moving segments ~2.3s each), then 56dps bars (shorter ~1.2s each).
    % We can split by identifying the duration change.

    fprintf('\n  --- Duration-based speed classification ---\n');
    if ~isempty(moving_durs)
        slow_count = sum(moving_durs > 1500);
        fast_count = sum(moving_durs <= 1500);
        fprintf('  Moving segments > 1.5s (28dps): %d\n', slow_count);
        fprintf('  Moving segments <= 1.5s (56dps): %d\n', fast_count);
    end

    %% 2f. Verify LUT directions
    fprintf('\n  --- LUT Direction Verification ---\n');
    slow_func_mask = (ce.func_order == 3 | ce.func_order == 4);
    n_slow = sum(slow_func_mask);
    fprintf('  Slow bar conditions (func 3 or 4): %d\n', n_slow);

    if n_slow == 16
        [lut_dirs, ~, lut_pats, lut_funs] = ...
            verify_lut_directions(Tbl, ce.pattern_order, ce.func_order, plot_order);

        fprintf('  LUT directions by data row:\n');
        fprintf('  %-8s %-8s %-8s %-12s\n', 'DataRow', 'Pattern', 'Func', 'Direction');
        for i = 1:16
            fprintf('  %-8d %-8d %-8d %-12.1f\n', i, lut_pats(i), lut_funs(i), lut_dirs(i));
        end

        % Check for pattern 9/10 (experiment patterns mapped to 0015/0016)
        fprintf('\n  --- Pattern 9/10 check (0015/0016 issue) ---\n');
        pat9_mask = lut_pats == 9;
        pat10_mask = lut_pats == 10;
        if any(pat9_mask)
            fprintf('  Pattern 9 rows: %s, directions: %s\n', ...
                mat2str(find(pat9_mask)'), mat2str(lut_dirs(pat9_mask)'));
        end
        if any(pat10_mask)
            fprintf('  Pattern 10 rows: %s, directions: %s\n', ...
                mat2str(find(pat10_mask)'), mat2str(lut_dirs(pat10_mask)'));
        end
    else
        fprintf('  WARNING: Only %d slow bar conditions found (expected 16)\n', n_slow);
        unique_funcs = unique(ce.func_order);
        for uf = 1:numel(unique_funcs)
            fprintf('    func %d: %d entries\n', unique_funcs(uf), ...
                sum(ce.func_order == unique_funcs(uf)));
        end
    end

    %% 2g. Full currentExp ordering for reference
    fprintf('\n  --- Full currentExp stimulus order (cycle 1) ---\n');
    n_per_cycle = numel(ce.pattern_order) / 3;
    if mod(numel(ce.pattern_order), 3) == 0
        fprintf('  Entries per cycle: %d (total: %d, 3 reps)\n', ...
            round(n_per_cycle), numel(ce.pattern_order));
        fprintf('  %-5s  %-10s  %-10s\n', 'Idx', 'Pattern', 'Function');
        for i = 1:min(round(n_per_cycle), 50)
            fprintf('  %-5d  %-10d  %-10d\n', i, ce.pattern_order(i), ce.func_order(i));
        end
    else
        fprintf('  WARNING: pattern_order length (%d) not divisible by 3\n', ...
            numel(ce.pattern_order));
    end
end

fprintf('\n========================================================================\n');
fprintf('  DIAGNOSTIC COMPLETE\n');
fprintf('========================================================================\n');
