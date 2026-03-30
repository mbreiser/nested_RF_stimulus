% COMPARE_PD_FLASH_IMPACT  Assess whether 4 PD methods change bar flash alignment.
%
%   For the 25 late-experiment cells (which have bar flash data), this script:
%     1. Loads batch_results.mat and computes 4 PD methods per cell
%     2. Converts each method's PD from aligned coordinates to stimulus coordinates
%     3. Snaps each to the nearest 22.5 deg grid direction
%     4. Builds a comparison table flagging cells where the grid direction differs
%     5. Saves the table to pd_flash_impact_table.txt
%     6. Regenerates individual cell polar-with-traces figures with 4 PD lines
%
%   OUTPUT FILES (in <data_root>/figure_previews/):
%     pd_flash_impact_table.txt       — per-cell PD comparison table
%     pd_4method_<date_str>.png       — individual cell polar figures with 4 PD lines
%
%   See also COMPUTE_PD_FOUR_METHODS, ADD_FOUR_PD_LINES, BATCH_ANALYZE_1DRF

%% Setup
data_root   = '/Users/reiserm/Documents/ttl_1DRF';
preview_dir = fullfile(data_root, 'figure_previews');
if ~isfolder(preview_dir), mkdir(preview_dir); end

% Load batch results
load(fullfile(data_root, 'population_results', 'batch_results.mat'), 'results');
n_cells = numel(results);
fprintf('Loaded %d cells from batch_results.mat\n', n_cells);

% Load LUT (same as batch pipeline)
script_dir = fileparts(mfilename('fullpath'));
lut_path = fullfile(fileparts(script_dir), 'src', 'analysis', 'protocol2', 'bar_lut.mat');
S_lut = load(lut_path, 'Tbl');
Tbl = S_lut.Tbl;

% Processing options (same as batch pipeline)
plot_order = [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16];
sweep_opts.baseline_range = [1000 9000];
sweep_opts.stim_trim_end  = 7000;
sweep_opts.percentile     = 98;
pattern_offset = 2;
on_threshold = 129;

%% ====================================================================
%%  PART 1: PD Comparison Table
%% ====================================================================

fprintf('\n=== Computing 4 PD methods for %d cells ===\n', n_cells);

% Preallocate table arrays
date_strs    = cell(n_cells, 1);
group_strs   = cell(n_cells, 1);
pd_stored    = NaN(n_cells, 1);       % stored VecNN PD (stimulus coords)
grid_pd      = NaN(n_cells, 4);       % snapped stimulus PD per method
stim_pd      = NaN(n_cells, 4);       % continuous stimulus PD per method
any_changed  = false(n_cells, 1);     % flag: any method differs from VecNN

method_names = {'VecNN', 'SymNN', 'VecInterp', 'SymInterp'};
method_fields = {'vec_nn', 'sym_nn', 'vec_interp', 'sym_interp'};

for k = 1:n_cells
    date_strs{k}  = results(k).date_str;
    group_strs{k} = results(k).group;
    pd_stored(k)  = results(k).pd_direction;

    % Get aligned tuning data (VecNN PD is at 90 deg in aligned coords)
    d = results(k).max_v_aligned;
    angles_16    = d(:, 1);
    responses_16 = d(:, 2);

    % Compute 4 PD methods
    methods = compute_pd_four_methods(angles_16, responses_16);

    % Convert each method's aligned PD to stimulus coordinates
    % In aligned coords: VecNN PD = 90 deg
    % Stimulus PD = aligned_pd + (stored_pd - 90)
    offset = pd_stored(k) - 90;

    for m = 1:4
        aligned_pd_deg = methods.(method_fields{m}).pd_deg;
        stim_pd(k, m) = mod(aligned_pd_deg + offset, 360);
        grid_pd(k, m) = mod(round(stim_pd(k, m) / 22.5) * 22.5, 360);
    end

    % Flag if any method's grid PD differs from VecNN grid PD
    any_changed(k) = any(grid_pd(k, 2:4) ~= grid_pd(k, 1));
end

%% Print and save table
fprintf('\n');
fprintf('==========================================================================\n');
fprintf('  PD Method Comparison — Bar Flash Impact (25 late cells)\n');
fprintf('==========================================================================\n');
header = sprintf('%-5s %-18s %-14s %-9s | %-9s %-9s %-9s %-9s | %-7s', ...
    '#', 'Date', 'Group', 'StoredPD', 'VecNN', 'SymNN', 'VecIntrp', 'SymIntrp', 'Changed');
fprintf('%s\n', header);
fprintf('%s\n', repmat('-', 1, numel(header)));

for k = 1:n_cells
    % Mark changed grid values with asterisk
    marks = repmat({' '}, 1, 4);
    for m = 2:4
        if grid_pd(k, m) ~= grid_pd(k, 1)
            marks{m} = '*';
        end
    end

    fprintf('%-5d %-18s %-14s %6.1f    | %6.1f%s  %6.1f%s  %6.1f%s  %6.1f%s  | %s\n', ...
        k, date_strs{k}, group_strs{k}, pd_stored(k), ...
        grid_pd(k,1), marks{1}, grid_pd(k,2), marks{2}, ...
        grid_pd(k,3), marks{3}, grid_pd(k,4), marks{4}, ...
        tern(any_changed(k), 'YES', ''));
end

fprintf('%s\n', repmat('-', 1, numel(header)));

%% Summary statistics
n_changed_per_method = zeros(1, 4);
for m = 1:4
    n_changed_per_method(m) = sum(grid_pd(:, m) ~= grid_pd(:, 1));
end

fprintf('\nSummary: Cells whose 22.5 deg grid PD differs from VecNN:\n');
for m = 2:4
    fprintf('  %s: %d / %d cells (%.0f%%)\n', method_names{m}, ...
        n_changed_per_method(m), n_cells, 100*n_changed_per_method(m)/n_cells);
end
fprintf('  Any method: %d / %d cells (%.0f%%)\n', ...
    sum(any_changed), n_cells, 100*sum(any_changed)/n_cells);

% Per-group breakdown
fprintf('\nPer-group breakdown (any method changes grid PD):\n');
for g = ["on_control", "on_ttl", "off_control", "off_ttl"]
    mask = strcmp(group_strs, g);
    n_g = sum(mask);
    n_c = sum(any_changed(mask));
    fprintf('  %-14s: %d / %d changed\n', g, n_c, n_g);
end

%% Save table to file
table_file = fullfile(preview_dir, 'pd_flash_impact_table.txt');
fid = fopen(table_file, 'w');
fprintf(fid, 'PD Method Comparison — Bar Flash Impact (%d late cells)\n', n_cells);
fprintf(fid, 'Generated: %s\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, 'Coordinate conversion:\n');
fprintf(fid, '  stim_pd = aligned_pd + (stored_pd - 90)\n');
fprintf(fid, '  grid_pd = round(stim_pd / 22.5) * 22.5  (mod 360)\n');
fprintf(fid, '  Changed = grid PD differs from VecNN grid PD\n');
fprintf(fid, '  * = this method''s grid PD differs from VecNN\n\n');

fprintf(fid, '%s\n', header);
fprintf(fid, '%s\n', repmat('-', 1, numel(header)));

for k = 1:n_cells
    marks = repmat({' '}, 1, 4);
    for m = 2:4
        if grid_pd(k, m) ~= grid_pd(k, 1)
            marks{m} = '*';
        end
    end
    fprintf(fid, '%-5d %-18s %-14s %6.1f    | %6.1f%s  %6.1f%s  %6.1f%s  %6.1f%s  | %s\n', ...
        k, date_strs{k}, group_strs{k}, pd_stored(k), ...
        grid_pd(k,1), marks{1}, grid_pd(k,2), marks{2}, ...
        grid_pd(k,3), marks{3}, grid_pd(k,4), marks{4}, ...
        tern(any_changed(k), 'YES', ''));
end

fprintf(fid, '%s\n\n', repmat('-', 1, numel(header)));

fprintf(fid, 'Summary: Cells whose 22.5 deg grid PD differs from VecNN:\n');
for m = 2:4
    fprintf(fid, '  %s: %d / %d cells (%.0f%%)\n', method_names{m}, ...
        n_changed_per_method(m), n_cells, 100*n_changed_per_method(m)/n_cells);
end
fprintf(fid, '  Any method: %d / %d cells (%.0f%%)\n', ...
    sum(any_changed), n_cells, 100*sum(any_changed)/n_cells);

fprintf(fid, '\nPer-group breakdown (any method changes grid PD):\n');
for g = ["on_control", "on_ttl", "off_control", "off_ttl"]
    mask = strcmp(group_strs, g);
    n_g = sum(mask);
    n_c = sum(any_changed(mask));
    fprintf(fid, '  %-14s: %d / %d changed\n', g, n_c, n_g);
end

fprintf(fid, '\nDetailed continuous PD (stimulus coords, before grid snap):\n');
fprintf(fid, '%-5s %-18s | %10s %10s %10s %10s\n', ...
    '#', 'Date', 'VecNN', 'SymNN', 'VecIntrp', 'SymIntrp');
fprintf(fid, '%s\n', repmat('-', 1, 70));
for k = 1:n_cells
    fprintf(fid, '%-5d %-18s | %10.1f %10.1f %10.1f %10.1f\n', ...
        k, date_strs{k}, stim_pd(k,1), stim_pd(k,2), stim_pd(k,3), stim_pd(k,4));
end

fclose(fid);
fprintf('\nTable saved to: %s\n', table_file);

%% ====================================================================
%%  PART 2: Individual Cell Polar Figures with 4 PD Lines
%% ====================================================================

fprintf('\n=== Generating individual cell polar figures with 4 PD lines ===\n');

for k = 1:n_cells
    date_str = results(k).date_str;
    exp_folder = fullfile(data_root, results(k).folder);
    fprintf('[%d/%d] %s (%s) ...', k, n_cells, date_str, results(k).group);

    try
        % Save/restore working directory (load_protocol2_data uses cd)
        orig_dir = pwd;

        % Load raw experiment data
        [~, ~, Log, ~, ~] = load_protocol2_data(exp_folder);

        f_data   = Log.ADC.Volts(1, :);
        v_data   = Log.ADC.Volts(2, :) * 10;
        median_v = median(v_data);

        ce = load(fullfile(exp_folder, 'currentExp.mat'), ...
            'pattern_order', 'func_order', 'metadata');
        metadata = ce.metadata;

        % Parse bar sweeps
        bar_data = parse_bar_data(f_data, v_data);

        [lut_directions, lut_orientations, lut_patterns, lut_functions] = ...
            verify_lut_directions(Tbl, ce.pattern_order, ce.func_order, plot_order);

        % Compute bar sweep responses
        max_v = compute_bar_sweep_responses(bar_data, plot_order, sweep_opts);

        % Determine ON/OFF label
        is_on = metadata.Frame > on_threshold;
        if is_on
            on_off = 'ON';
        else
            on_off = 'OFF';
        end

        % Create polar figure (ring of traces + central polar)
        change_str = '';
        if any_changed(k), change_str = ' [PD CHANGED]'; end
        polar_title = sprintf('28 dps - %s - %s - %s%s', ...
            strrep(date_str, '_', '-'), ...
            strrep(metadata.Strain, '_', ' '), ...
            on_off, change_str);

        fig = plot_slow_bar_sweep_polar(bar_data, max_v, lut_directions, ...
            plot_order, median_v, polar_title);

        % Add 4 PD direction lines
        % pd_stim_deg_4 = [VecNN, SymNN, VecInterp, SymInterp] in stimulus degrees
        pd_stim_deg_4 = stim_pd(k, :);

        % Find max rho for line extent
        max_rho = max(max_v);

        add_four_pd_lines(fig, pd_stim_deg_4, max_rho);

        % Save PNG (use folder name to avoid duplicate date collisions)
        png_path = fullfile(preview_dir, sprintf('pd_4method_%s.png', results(k).folder));
        exportgraphics(fig, png_path, 'Resolution', 200);
        close(fig);

        cd(orig_dir);

        fprintf(' saved\n');

    catch ME
        cd(orig_dir);
        fprintf(' ERROR: %s\n', ME.message);
        continue;
    end
end

fprintf('\n=== Done. %d PNGs saved to %s ===\n', n_cells, preview_dir);


%% ========================= Local Helpers ============================

function s = tern(cond, a, b)
% TERN  Inline ternary operator.
    if cond
        s = a;
    else
        s = b;
    end
end
