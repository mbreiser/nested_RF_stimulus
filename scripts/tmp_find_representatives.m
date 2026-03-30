% Temporary script to find representative cells closest to group mean DSI
addpath(genpath('src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

late = load('/Users/reiserm/Documents/ttl_1DRF/population_results/batch_results.mat');
early = load('/Users/reiserm/Documents/ttl_1DRF/pre-bar-flash/population_results/batch_results_pre_bf.mat');

groups = {'ON_ctrl', 'ON_ttl', 'OFF_ctrl', 'OFF_ttl'};
for g = 1:4
    grp = groups{g};
    dsi_vals = [];
    cell_ids = {};
    sources = {};

    for i = 1:numel(late.results)
        r = late.results(i);
        match = false;
        if r.is_on && startsWith(grp, 'ON') && ~r.is_ttl && endsWith(grp, 'ctrl'), match = true; end
        if r.is_on && startsWith(grp, 'ON') && r.is_ttl && endsWith(grp, 'ttl'), match = true; end
        if ~r.is_on && startsWith(grp, 'OFF') && ~r.is_ttl && endsWith(grp, 'ctrl'), match = true; end
        if ~r.is_on && startsWith(grp, 'OFF') && r.is_ttl && endsWith(grp, 'ttl'), match = true; end
        if match
            dsi_vals(end+1) = r.dsi_vector;
            cell_ids{end+1} = r.date_str;
            sources{end+1} = 'late';
        end
    end

    for i = 1:numel(early.results)
        r = early.results(i);
        match = false;
        if r.is_on && startsWith(grp, 'ON') && ~r.is_ttl && endsWith(grp, 'ctrl'), match = true; end
        if r.is_on && startsWith(grp, 'ON') && r.is_ttl && endsWith(grp, 'ttl'), match = true; end
        if ~r.is_on && startsWith(grp, 'OFF') && ~r.is_ttl && endsWith(grp, 'ctrl'), match = true; end
        if ~r.is_on && startsWith(grp, 'OFF') && r.is_ttl && endsWith(grp, 'ttl'), match = true; end
        if match
            dsi_vals(end+1) = r.dsi_vector;
            cell_ids{end+1} = r.date_str;
            sources{end+1} = 'early';
        end
    end

    mean_dsi = mean(dsi_vals);
    [~, closest_idx] = min(abs(dsi_vals - mean_dsi));

    fprintf('\n%s (n=%d): mean DSI = %.3f\n', grp, numel(dsi_vals), mean_dsi);
    fprintf('  BEST: %s (source=%s, DSI=%.3f, delta=%.4f)\n', ...
        cell_ids{closest_idx}, sources{closest_idx}, dsi_vals(closest_idx), ...
        abs(dsi_vals(closest_idx) - mean_dsi));

    [sorted_delta, sort_idx] = sort(abs(dsi_vals - mean_dsi));
    for j = 1:min(5, numel(sort_idx))
        k = sort_idx(j);
        fprintf('    %d. %s (%s) DSI=%.3f delta=%.4f\n', j, cell_ids{k}, sources{k}, dsi_vals(k), sorted_delta(j));
    end
end
