function [X, t_ms] = step_snips(v, k_on, win_ms, fs)
%STEP_SNIPS  Snippets of a 10 kHz trace around onsets, at 1 kHz (10-sample block means, as pharm_protocol.decimate).
%
%   [X, T_MS] = STEP_SNIPS(V, K_ON, WIN_MS) returns X (numel(K_ON) x T) with the mean of V over each
%   1-ms bin [t, t+1) ms from WIN_MS(1) to WIN_MS(2) after each onset sample K_ON, and T_MS the bin
%   starts. Rows of onsets that are NaN or too close to the ends are NaN. FS default 10000.
%   28 October 2026.
if nargin < 4, fs = 10000; end
spb = fs / 1000; t_ms = win_ms(1):win_ms(2) - 1; nb = numel(t_ms);
X = nan(numel(k_on), nb);
off = (t_ms(1) * spb) + (0:nb * spb - 1);
for i = 1:numel(k_on)
    if isnan(k_on(i)), continue; end
    idx = k_on(i) + off;
    if idx(1) < 1 || idx(end) > numel(v), continue; end
    X(i, :) = mean(reshape(v(idx), spb, nb), 1);
end
end
