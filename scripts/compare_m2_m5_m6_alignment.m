%% compare_m2_m5_m6_alignment.m
%  Compare three RF-center methods on both PD and ortho axes for all 25 cells.
%
%  M2 — argmax (position of peak amplitude)
%  M5 — FWHM bump centroid (center-of-mass within half-max region)
%  M6 — 68%-area centroid (center-of-mass within smallest contiguous
%        region around peak containing >=68% of total area)
%
%  Outputs a formatted table and summary statistics.

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src')));

%% Load data
data_root = '/Users/reiserm/Documents/ttl_1DRF';
load(fullfile(data_root, 'population_results', 'batch_results.mat'), 'results');
n_cells = numel(results);
fprintf('Loaded %d cells\n\n', n_cells);

%% Compute M2, M5, M6 for each cell on both axes

% Preallocate
pd_m2  = nan(n_cells, 1);   pd_m5  = nan(n_cells, 1);   pd_m6  = nan(n_cells, 1);
pd_m5f = nan(n_cells, 1);   pd_m6f = nan(n_cells, 1);  % fractional
od_m2  = nan(n_cells, 1);   od_m5  = nan(n_cells, 1);   od_m6  = nan(n_cells, 1);
od_m5f = nan(n_cells, 1);   od_m6f = nan(n_cells, 1);
pd_bw5 = nan(n_cells, 1);   pd_bw6 = nan(n_cells, 1);   % bump widths
od_bw5 = nan(n_cells, 1);   od_bw6 = nan(n_cells, 1);

for i = 1:n_cells
    % --- PD axis ---
    A_pd = max(results(i).peak_amplitudes, 0);
    [~, pd_m2(i)] = max(A_pd);

    m5_pd = compute_m5_centroid(A_pd);
    pd_m5f(i) = m5_pd.centroid;
    pd_m5(i)  = m5_pd.centroid_int;
    pd_bw5(i) = m5_pd.bump_width;

    m6_pd = compute_m6_centroid(A_pd, 0.68);
    pd_m6f(i) = m6_pd.centroid;
    pd_m6(i)  = m6_pd.centroid_int;
    pd_bw6(i) = m6_pd.bump_width;

    % --- Ortho axis ---
    A_od = max(results(i).ortho_peak_amplitudes, 0);
    [~, od_m2(i)] = max(A_od);

    m5_od = compute_m5_centroid(A_od);
    od_m5f(i) = m5_od.centroid;
    od_m5(i)  = m5_od.centroid_int;
    od_bw5(i) = m5_od.bump_width;

    m6_od = compute_m6_centroid(A_od, 0.68);
    od_m6f(i) = m6_od.centroid;
    od_m6(i)  = m6_od.centroid_int;
    od_bw6(i) = m6_od.bump_width;
end

%% Build comparison table

% Group labels
groups = cell(n_cells, 1);
for i = 1:n_cells
    if results(i).is_on && ~results(i).is_ttl
        groups{i} = 'ON ctrl';
    elseif results(i).is_on && results(i).is_ttl
        groups{i} = 'ON TTL';
    elseif ~results(i).is_on && ~results(i).is_ttl
        groups{i} = 'OFF ctrl';
    else
        groups{i} = 'OFF TTL';
    end
end

% Print table header
fprintf('=== PD Axis: M2 vs M5 vs M6 per cell ===\n\n');
fprintf('%-4s %-10s %-20s  M2  M5(frac)  M5  M6(frac)  M6  BW5  BW6  M2==M5  M2==M6  M5==M6\n', ...
    '#', 'Group', 'Experiment');
fprintf('%s\n', repmat('-', 1, 110));

n_m2_m5_diff_pd = 0;
n_m2_m6_diff_pd = 0;
n_m5_m6_diff_pd = 0;

for i = 1:n_cells
    m2m5 = '  yes ';
    if pd_m2(i) ~= pd_m5(i), m2m5 = '  NO  '; n_m2_m5_diff_pd = n_m2_m5_diff_pd + 1; end
    m2m6 = '  yes ';
    if pd_m2(i) ~= pd_m6(i), m2m6 = '  NO  '; n_m2_m6_diff_pd = n_m2_m6_diff_pd + 1; end
    m5m6 = '  yes ';
    if pd_m5(i) ~= pd_m6(i), m5m6 = '  NO  '; n_m5_m6_diff_pd = n_m5_m6_diff_pd + 1; end

    fprintf('%-4d %-10s %-20s  %2d   %5.2f    %2d   %5.2f    %2d   %2d   %2d  %s  %s  %s\n', ...
        i, groups{i}, results(i).date_str, ...
        pd_m2(i), pd_m5f(i), pd_m5(i), pd_m6f(i), pd_m6(i), ...
        pd_bw5(i), pd_bw6(i), m2m5, m2m6, m5m6);
end

fprintf('\nPD summary: M2≠M5 in %d/%d cells, M2≠M6 in %d/%d cells, M5≠M6 in %d/%d cells\n', ...
    n_m2_m5_diff_pd, n_cells, n_m2_m6_diff_pd, n_cells, n_m5_m6_diff_pd, n_cells);

% Ortho table
fprintf('\n\n=== Ortho Axis: M2 vs M5 vs M6 per cell ===\n\n');
fprintf('%-4s %-10s %-20s  M2  M5(frac)  M5  M6(frac)  M6  BW5  BW6  M2==M5  M2==M6  M5==M6\n', ...
    '#', 'Group', 'Experiment');
fprintf('%s\n', repmat('-', 1, 110));

n_m2_m5_diff_od = 0;
n_m2_m6_diff_od = 0;
n_m5_m6_diff_od = 0;

for i = 1:n_cells
    m2m5 = '  yes ';
    if od_m2(i) ~= od_m5(i), m2m5 = '  NO  '; n_m2_m5_diff_od = n_m2_m5_diff_od + 1; end
    m2m6 = '  yes ';
    if od_m2(i) ~= od_m6(i), m2m6 = '  NO  '; n_m2_m6_diff_od = n_m2_m6_diff_od + 1; end
    m5m6 = '  yes ';
    if od_m5(i) ~= od_m6(i), m5m6 = '  NO  '; n_m5_m6_diff_od = n_m5_m6_diff_od + 1; end

    fprintf('%-4d %-10s %-20s  %2d   %5.2f    %2d   %5.2f    %2d   %2d   %2d  %s  %s  %s\n', ...
        i, groups{i}, results(i).date_str, ...
        od_m2(i), od_m5f(i), od_m5(i), od_m6f(i), od_m6(i), ...
        od_bw5(i), od_bw6(i), m2m5, m2m6, m5m6);
end

fprintf('\nOrtho summary: M2≠M5 in %d/%d cells, M2≠M6 in %d/%d cells, M5≠M6 in %d/%d cells\n', ...
    n_m2_m5_diff_od, n_cells, n_m2_m6_diff_od, n_cells, n_m5_m6_diff_od, n_cells);

%% Cross-axis summary: how far each method moves from M2
fprintf('\n\n=== Summary: mean |shift| from M2 (in positions) ===\n\n');
fprintf('Axis    Method   MeanShift   MaxShift   Cells_shifted\n');
fprintf('------  ------   ---------   --------   -------------\n');
fprintf('PD      M5       %6.2f       %d          %d/%d\n', ...
    mean(abs(pd_m5 - pd_m2)), max(abs(pd_m5 - pd_m2)), n_m2_m5_diff_pd, n_cells);
fprintf('PD      M6       %6.2f       %d          %d/%d\n', ...
    mean(abs(pd_m6 - pd_m2)), max(abs(pd_m6 - pd_m2)), n_m2_m6_diff_pd, n_cells);
fprintf('Ortho   M5       %6.2f       %d          %d/%d\n', ...
    mean(abs(od_m5 - od_m2)), max(abs(od_m5 - od_m2)), n_m2_m5_diff_od, n_cells);
fprintf('Ortho   M6       %6.2f       %d          %d/%d\n', ...
    mean(abs(od_m6 - od_m2)), max(abs(od_m6 - od_m2)), n_m2_m6_diff_od, n_cells);

%% Bump width comparison
fprintf('\n\n=== Bump widths by method and group ===\n\n');
fprintf('Group       n   PD_BW5  PD_BW6  OD_BW5  OD_BW6\n');
fprintf('----------  --  ------  ------  ------  ------\n');

group_names = {'ON ctrl', 'ON TTL', 'OFF ctrl', 'OFF TTL'};
for g = 1:4
    mask = strcmp(groups, group_names{g});
    fprintf('%-10s  %2d  %5.1f   %5.1f   %5.1f   %5.1f\n', ...
        group_names{g}, sum(mask), ...
        mean(pd_bw5(mask)), mean(pd_bw6(mask)), ...
        mean(od_bw5(mask)), mean(od_bw6(mask)));
end
fprintf('%-10s  %2d  %5.1f   %5.1f   %5.1f   %5.1f\n', ...
    'ALL', n_cells, mean(pd_bw5), mean(pd_bw6), mean(od_bw5), mean(od_bw6));


%% ====================================================================
%  Local functions
%  ====================================================================

function m5 = compute_m5_centroid(A)
% COMPUTE_M5_CENTROID  FWHM bump centroid.
%   Walk from peak in both directions while amp >= peak/2.
%   Center-of-mass within that region.

    n_pos = numel(A);
    [peak_val, peak_pos] = max(A);

    if peak_val <= 0
        m5.centroid     = 6;
        m5.centroid_int = 6;
        m5.bump_range   = [1, n_pos];
        m5.bump_width   = n_pos;
        return;
    end

    threshold = peak_val / 2;

    left = peak_pos;
    while left > 1 && A(left - 1) >= threshold
        left = left - 1;
    end
    right = peak_pos;
    while right < n_pos && A(right + 1) >= threshold
        right = right + 1;
    end

    bump = left:right;
    A_bump = A(bump);
    m5.centroid     = sum(A_bump .* bump) / sum(A_bump);
    m5.centroid_int = round(m5.centroid);
    m5.centroid_int = max(1, min(n_pos, m5.centroid_int));
    m5.bump_range   = [left, right];
    m5.bump_width   = right - left + 1;
end


function m6 = compute_m6_centroid(A, area_fraction)
% COMPUTE_M6_CENTROID  Area-fraction bump centroid.
%   Grow a contiguous region outward from the peak until it contains
%   >= area_fraction (e.g. 0.68) of the total positive area.
%   Center-of-mass within that region.
%
%   Algorithm:
%     1. Start with the peak position only.
%     2. At each step, look at the two neighbors just outside the current
%        region (left-1 and right+1).  Add whichever has the larger
%        amplitude (greedy expansion toward mass).
%     3. Stop when cumulative area >= area_fraction * total_area,
%        or the region spans all positions.
%   This produces a narrower region than FWHM for broad bumps, staying
%   closer to the peak.

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
        % Look at candidates
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
