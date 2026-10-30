function data = parse_bar_data_pre_bf(f_data, v_data, speed)
% PARSE_BAR_DATA_PRE_BF  Parse bar sweep data from pre-bar-flash protocol.
%
%   DATA = PARSE_BAR_DATA_PRE_BF(F_DATA, V_DATA) extracts slow (28 dps) bar
%   sweep voltage traces from the 2-speed pre-bar-flash protocol (summer 2025).
%
%   DATA = PARSE_BAR_DATA_PRE_BF(F_DATA, V_DATA, SPEED) extracts bars at the
%   specified speed: 'slow' (28 dps, default) or 'medium' (56 dps).
%
%   This parser handles both pre-Oct-15 experiments (no intra-cycle grey gaps)
%   and post-Oct-15 experiments (with 3-4s intra-cycle gaps, extra speeds, and
%   bar flashes). Extracts only bars matching the requested speed.
%
%   Protocol structure (per cycle, pre-Oct-15):
%     10s grey → 4px flashes → 6px flashes → 28 dps bars → 56 dps bars
%     (bars run continuously — no grey gaps between speed blocks)
%
%   Protocol structure (per cycle, post-Oct-15):
%     10s grey → flashes → [3s gap] → 28 dps → [3s gap] → 56 dps →
%     [3s gap] → 168 dps → [3s gap] → bar flashes → [3s gap] → bar flashes
%
%   INPUTS:
%     f_data - Frame signal (1xM), 10 kHz. 0 = grey screen, >0 = pattern ID
%     v_data - Voltage signal (1xM), 10 kHz, in mV
%     speed  - (Optional) 'slow' for 28 dps (default), 'medium' for 56 dps
%
%   OUTPUT:
%     data - 16x4 cell array. Rows = bar directions in presentation order.
%            Columns 1-3 = individual cycle traces, column 4 = mean.
%            Each trace includes 9000 samples pre- and post-stimulus padding.
%
%   ALGORITHM:
%     1. Detect cycle boundaries from >= 8s grey-screen gaps
%     2. Within each cycle, find all frame transitions (|diff| > 9)
%     3. Build segments between transitions
%     4. Classify segments: moving (frame std > 5) AND duration in speed range
%     5. Extract 16 matching segments per cycle with 9000-sample padding
%     6. Average across cycles into column 4
%
%   NOTES:
%     Duration-based speed classification:
%       28 dps bars: ~21000 samples (2.1s) → 'slow' (dur > 15000)
%       56 dps bars: ~10700 samples (1.07s) → 'medium' (8000 < dur < 15000)
%       168 dps bars: ~3600 samples (0.36s) → not supported
%       Flash segments: <1000 samples → excluded
%
%   See also PARSE_BAR_DATA, COMPUTE_BAR_SWEEP_RESPONSES, BATCH_ANALYZE_PRE_BAR_FLASH

    %% Parameters
    if nargin < 3 || isempty(speed), speed = 'slow'; end

    GREY_GAP_MIN   = 80000;  % 8s at 10 kHz — separates cycle-boundary gaps from intra-cycle gaps
    TRANS_THR      = 9;      % Min frame value change to detect segment boundary
    MOVING_STD_THR = 5;      % Frame std threshold: moving bar > 5, static bar ~ 0
    SLOW_DUR_MIN   = 15000;  % 1.5s — cleanly separates 28 dps (>2s) from 56 dps (<1.1s)
    MED_DUR_MIN    = 8000;   % 0.8s — lower bound for 56 dps bars (~10700 samples)
    MED_DUR_MAX    = 15000;  % 1.5s — upper bound for 56 dps bars (excludes 28 dps)
    PADDING        = 9000;   % Samples before/after bar onset (matches parse_bar_data)
    N_DIRS         = 16;     % Expected bar directions per cycle

    %% 1. Find cycle boundaries from long grey-screen gaps
    zero_mask = f_data == 0;
    d = diff([0, zero_mask, 0]);
    gap_starts = find(d == 1);
    gap_ends   = find(d == -1) - 1;
    gap_lens   = gap_ends - gap_starts + 1;

    long_mask = gap_lens >= GREY_GAP_MIN;
    long_idx  = find(long_mask);
    n_long    = numel(long_idx);

    if n_long < 2
        error('parse_bar_data_pre_bf:noCycles', ...
            'Found %d grey gaps >= 8s — need at least 2 to define cycles.', n_long);
    end

    n_cycles = n_long - 1;

    %% 2. For each cycle, find slow moving bar segments
    all_slow_segs = cell(1, n_cycles);  % each: Kx2 [abs_start, abs_end]

    for c = 1:n_cycles
        % Cycle range: from end of gap c to start of gap c+1
        c_start = gap_ends(long_idx(c)) + 1;
        c_end   = gap_starts(long_idx(c + 1)) - 1;

        f_cycle = f_data(c_start:c_end);

        % Find all frame transitions within this cycle
        df = abs(diff(double(f_cycle)));
        trans_local = find(df > TRANS_THR);

        if numel(trans_local) < 2
            warning('parse_bar_data_pre_bf:fewTransitions', ...
                'Cycle %d: only %d transitions found.', c, numel(trans_local));
            all_slow_segs{c} = zeros(0, 2);
            continue;
        end

        % Convert to absolute indices.
        % trans_local(s) means f_cycle(s) -> f_cycle(s+1) has a big change.
        % In absolute coords: f_data(c_start + s - 1) -> f_data(c_start + s).
        % The segment AFTER transition s starts at c_start + trans_local(s).
        trans_abs = c_start + trans_local;

        % Build segments between consecutive transitions and classify.
        % Include a final segment from last transition to cycle end, since
        % the last bar in a speed block may extend to the cycle boundary.
        n_trans = numel(trans_abs);
        seg_bounds = zeros(n_trans, 2);
        for s = 1:n_trans - 1
            seg_bounds(s, :) = [trans_abs(s), trans_abs(s + 1) - 1];
        end
        seg_bounds(n_trans, :) = [trans_abs(n_trans), c_end];  % last-to-boundary

        n_slow = 0;
        slow_list = zeros(N_DIRS + 4, 2);  % slight over-allocation for safety

        for s = 1:n_trans
            seg_start = seg_bounds(s, 1);
            seg_end   = seg_bounds(s, 2);
            seg_dur   = seg_end - seg_start + 1;

            if seg_dur < 1, continue; end

            seg_std = std(double(f_data(seg_start:seg_end)));

            if strcmp(speed, 'slow')
                is_target = seg_std > MOVING_STD_THR && seg_dur > SLOW_DUR_MIN;
            else  % 'medium'
                is_target = seg_std > MOVING_STD_THR && ...
                    seg_dur > MED_DUR_MIN && seg_dur < MED_DUR_MAX;
            end
            if is_target
                n_slow = n_slow + 1;
                if n_slow <= size(slow_list, 1)
                    slow_list(n_slow, :) = [seg_start, seg_end];
                end
            end
        end

        if n_slow ~= N_DIRS
            warning('parse_bar_data_pre_bf:segCount', ...
                'Cycle %d: found %d slow moving bars (expected %d).', ...
                c, n_slow, N_DIRS);
        end

        all_slow_segs{c} = slow_list(1:min(n_slow, size(slow_list, 1)), :);
    end

    %% 3. Validate: all cycles should have the same number of slow bars
    n_per_cycle = cellfun(@(x) size(x, 1), all_slow_segs);
    n_bars = min(n_per_cycle);

    if n_bars == 0
        error('parse_bar_data_pre_bf:noBars', ...
            'No slow moving bars found in any cycle.');
    end
    if n_bars < N_DIRS
        warning('parse_bar_data_pre_bf:barCount', ...
            'Using %d bars per cycle (min across cycles; expected %d).', ...
            n_bars, N_DIRS);
    end

    %% 4. Extract voltage traces with pre/post padding
    M = numel(v_data);
    data = cell(n_bars, n_cycles + 1);

    for c = 1:n_cycles
        segs = all_slow_segs{c};
        for b = 1:n_bars
            t_start = max(1, segs(b, 1) - PADDING);
            t_end   = min(M, segs(b, 2) + PADDING);
            data{b, c} = v_data(t_start:t_end);
        end
    end

    %% 5. Average across cycles into final column (matches parse_bar_data contract)
    for b = 1:n_bars
        % Collect per-cycle traces as column vectors
        d_cells = cellfun(@(x) x(:), data(b, 1:n_cycles), 'UniformOutput', false);
        lens = cellfun(@numel, d_cells);
        min_len = min(lens);

        % Trim to common length and average
        d_trimmed = cellfun(@(x) x(1:min_len), d_cells, 'UniformOutput', false);
        d_matrix  = horzcat(d_trimmed{:});
        data{b, n_cycles + 1} = mean(d_matrix, 2, 'omitnan')';
    end

end
