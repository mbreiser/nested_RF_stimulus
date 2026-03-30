function metrics = extract_temporal_metrics(traces_11xN, peak_pos, fs)
% EXTRACT_TEMPORAL_METRICS  Temporal analysis per position.
%
%   METRICS = EXTRACT_TEMPORAL_METRICS(TRACES_11XN, PEAK_POS, FS)
%   extracts timing metrics from baseline-subtracted flash traces.
%
%   For each position, measures:
%     - rise_start: time when trace first reaches 10% of max response (ms)
%     - rise_time:  time from 10% to 50% of max (ms)
%     - decay_time: time from 80% to 20% of max (ms)
%     - time_to_90: time from flash onset to 90% of max (ms)
%
%   All times are relative to flash onset (sample 5001 at 10 kHz).
%   Delta variants report each metric minus the value at peak_pos.
%
%   INPUTS:
%     traces_11xN - 11 x N baseline-subtracted mean flash traces
%                   (rows = positions ordered ND to PD, columns = samples)
%     peak_pos    - Integer position of peak response (1-11)
%     fs          - Sampling rate in Hz (default: 10000)
%
%   OUTPUT:
%     metrics - Structure with fields:
%       .rise_start       - 1x11 rise start times (ms from flash onset)
%       .rise_time        - 1x11 rise times (10% to 50%, ms)
%       .decay_time       - 1x11 decay times (80% to 20%, ms)
%       .time_to_90       - 1x11 time from flash onset to 90% of peak (ms)
%       .rise_start_delta - 1x11 rise_start minus rise_start(peak_pos)
%       .rise_time_delta  - 1x11 rise_time minus rise_time(peak_pos)
%       .decay_time_delta - 1x11 decay_time minus decay_time(peak_pos)
%       .time_to_90_delta - 1x11 time_to_90 minus time_to_90(peak_pos)
%       .peak_response    - 1x11 max response at each position (mV)
%
%   REFERENCE:
%     Gruntman et al. (eLife 2019): rise start (10% of max), rise time
%     (10% to 50%), decay time (80% to 20%). Reported as difference from
%     the central (peak) position.
%
%   See also BATCH_ANALYZE_1DRF, COMPUTE_RF_CENTROID

    if nargin < 3, fs = 10000; end

    n_pos = size(traces_11xN, 1);
    onset_sample = 5001;  % flash onset

    % Preallocate with NaN
    rise_start    = NaN(1, n_pos);
    rise_time     = NaN(1, n_pos);
    decay_time    = NaN(1, n_pos);
    time_to_90    = NaN(1, n_pos);
    peak_response = NaN(1, n_pos);

    for pos = 1:n_pos
        trace = traces_11xN(pos, :);

        % Skip all-NaN rows (from alignment padding)
        if all(isnan(trace))
            continue;
        end

        % Response portion: from flash onset to end of trace
        if onset_sample > numel(trace)
            continue;
        end
        resp = trace(onset_sample:end);
        peak_val = max(resp);
        peak_response(pos) = peak_val;

        if peak_val <= 0
            continue;
        end

        % Thresholds
        thresh_10 = 0.10 * peak_val;
        thresh_50 = 0.50 * peak_val;
        thresh_80 = 0.80 * peak_val;
        thresh_90 = 0.90 * peak_val;
        thresh_20 = 0.20 * peak_val;

        % Rise start: first sample >= 10% of peak
        idx_10 = find(resp >= thresh_10, 1, 'first');
        if isempty(idx_10)
            continue;
        end
        rise_start(pos) = (idx_10 - 1) / fs * 1000;  % ms from onset

        % Rise time: 10% to 50%
        idx_50 = find(resp >= thresh_50, 1, 'first');
        if ~isempty(idx_50) && idx_50 >= idx_10
            rise_time(pos) = (idx_50 - idx_10) / fs * 1000;
        end

        % Time to 90%: flash onset to 90% of peak
        idx_90 = find(resp >= thresh_90, 1, 'first');
        if ~isempty(idx_90)
            time_to_90(pos) = (idx_90 - 1) / fs * 1000;  % ms from onset
        end

        % Decay time: 80% to 20% after the peak sample
        [~, idx_peak] = max(resp);
        if idx_peak < numel(resp)
            decay_trace = resp(idx_peak:end);
            idx_80_decay = find(decay_trace <= thresh_80, 1, 'first');
            idx_20_decay = find(decay_trace <= thresh_20, 1, 'first');
            if ~isempty(idx_80_decay) && ~isempty(idx_20_decay) && idx_20_decay > idx_80_decay
                decay_time(pos) = (idx_20_decay - idx_80_decay) / fs * 1000;
            end
        end
    end

    % Pack raw metrics
    metrics.rise_start    = rise_start;
    metrics.rise_time     = rise_time;
    metrics.decay_time    = decay_time;
    metrics.time_to_90    = time_to_90;
    metrics.peak_response = peak_response;

    % Delta metrics (relative to peak position)
    if peak_pos >= 1 && peak_pos <= n_pos
        metrics.rise_start_delta = rise_start - rise_start(peak_pos);
        metrics.rise_time_delta  = rise_time  - rise_time(peak_pos);
        metrics.decay_time_delta = decay_time - decay_time(peak_pos);
        metrics.time_to_90_delta = time_to_90 - time_to_90(peak_pos);
    else
        metrics.rise_start_delta = NaN(1, n_pos);
        metrics.rise_time_delta  = NaN(1, n_pos);
        metrics.decay_time_delta = NaN(1, n_pos);
        metrics.time_to_90_delta = NaN(1, n_pos);
    end

end
