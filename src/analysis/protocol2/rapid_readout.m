function R = rapid_readout(session_folder, varargin)
% RAPID_READOUT  PD, RF center and the Prot3b / Protocol 4 placement tables from a Protocol 2 session, in MATLAB, in seconds.
%
%   R = RAPID_READOUT(SESSION_FOLDER) reads the raw G4 TDMS log of a nested-RF
%   (Protocol 2) session with FAST_TDMS_FOLDER2STRUCT (no conversion to .mat
%   needed) and computes, as tools/rapid_prot2_locations.py of
%   mbreiser/ephys-standard-analysis does (revision of 26 Oct 2026, 15:05):
%     - the TUTL metrics: the snapped PD (vector sum of the 28 deg/s sweep
%       amplitudes, 98th percentile of the repeat-mean trace over motion +
%       150 ms minus the pre-motion baseline, over the 16 directions, snapped to
%       the nearest LUT direction; VECTOR_SUM_POLAR / FIND_PD_FROM_LUT
%       convention) and the M6 center bars (COMPUTE_M6_CENTROID, 68 %-area
%       centroid of the 99.5th-percentile 80 ms bar-flash amplitudes) of the
%       bar set perpendicular to the snapped PD (PD axis) and of the set
%       parallel to it (orthogonal axis; + = 90 deg clockwise from PD);
%     - the placement PD: the 28 and 56 deg/s vector PDs separately; their
%       circular mean when they differ by <= 'AgreeDeg' (pd_rule 'mean'), else
%       the 28 deg/s vector alone with a warning ('28_only'); a disagreement
%       above 'FlagDeg' is flagged (pd_flag_gt5). The window, the sites, the lattice and
%       the two bar sets of the placement frame use this PD, and the M6 center
%       on each axis is recomputed from those sets;
%   prints them with the margin of the 28 deg/s vector to the nearest LUT
%   boundary, draws the placement figure (PLOT_RAPID_READOUT) and writes the
%   tables of the Python read-out with the same file names and columns, so
%   that the two can be diffed (tools/compare_readout_outputs.py there):
%     <cell>_sites.csv, <cell>_lattice.csv, <cell>_prot3b_frames.csv,
%     <cell>_prot3b_blocks.csv, <cell>_prot4_frames.csv,
%     <cell>_prot4_sequence.csv, <cell>_locations.json, <cell>_locations.png
%   On a stage-1 session (sweeps only, no bar flashes: the two-stage Protocol
%   2) it prints the placement PD and writes <cell>_stage1.json,
%   <cell>_stage2_bars.csv (the 22 bars of the two placement-frame sets as
%   frames of the 88-frame bar pattern) and <cell>_stage2_sequence.csv.
%   Specifications: docs/briefs/JY_PROTOCOL_HANDOVER_20261026.md (that
%   repository). The PD-frame offset s is positive toward the PD (the
%   trailing / proximal side); v is the orthogonal offset; arena positions are
%   1-based (column, row) with row 1 at the top; frame 1 of both patterns is
%   the background.
%
%   Options (name/value):
%     'Cell'     label used in the printout and the file names (default: the
%                session folder name)
%     'Out'      output folder (default: <session_folder>/rapid_readout;
%                refused under /Volumes/reiserlab, the read-only archive)
%     'Seed'     seed of the random block and sequence orders (default:
%                CRC-32 of the session id, as the Python script; the
%                permutations themselves differ from Python's because the
%                generators differ, see the note)
%     'Metrics'  metrics.csv of a dataset runner: snapped PD, 28 deg/s vector
%                PD, TUTL M6 bars and bar sets are compared with the row of
%                this cell / session
%     'Figure'   draw the figure (default true); 'SavePng' write it (true)
%     'Stage2'   two-stage Protocol 2, call 2: the folder of the stage-2 log
%                (e.g. <session>/p2b or <session>/p2b/Log Files). The sweeps
%                come from SESSION_FOLDER (the stage-1 package), the bar flashes
%                from this log; its Patterns/ and Functions/ are taken from the
%                stage-2 folder or its parent when present, else from the
%                stage-1 package (the bar pattern must be the same file)
%     'P1Peak'   Protocol 1 ON peak, 1-based [col row] (from RAPID_P1_PEAK):
%                written to other_centers.protocol1_ON_peak with its distance
%                to the center and drawn in the figure; 'P1Summary' takes it
%                from a p1_summary.csv row of this cell instead
%     'Verbose'  print the per-file read times (false)
%     'Stage1Test' emulate a stage-1 session from a complete one by dropping
%                the bar-flash trials (false; the Python --stage1-test)
%     'Stage2Test' emulate a stage-2 session from a complete one: keep only the
%                bar flashes of the 22 bars of the two placement-frame sets
%                (false). The center must equal the full-session center.
%     'AgreeDeg' the 28 and 56 deg/s vectors are averaged when they differ by
%                at most this (10, as the Python script at commit c372a6e; the
%                26 Oct decision text says 5); 'FlagDeg' flags a disagreement
%                above this without changing the rule (5)
%
%   R holds everything (session structure, amplitudes, PD, placement, bar
%   sets, geometry, sites, lattice, tables, validation, timing) for the
%   session sequencer.
%
%   Segmentation note. Trials are taken from the G4 command log (Set Pattern
%   ID / Set Pattern Function ID / Start-Display) and the events from the
%   frame channel inside each trial, as pharm_protocol.py does, so the same
%   code serves the 2025 three-speed sessions and Jin Yong's five-speed
%   sessions. PARSE_BAR_DATA / PARSE_BAR_FLASH_DATA index the 3-s gray gaps
%   assuming the 2025 block order and do not segment the five-speed protocol.
%
%   Validation (docs/notes/MATLAB_FAST_READOUT_20261026.md): identical to the
%   Python read-out on 25 TUTL sessions (historic layout) and on five October
%   2026 sessions of Jin Yong's rig (ADC1/ADC2 layout), where the TUTL fields
%   also PASS against datasets/20261022_jy_rnai_rf/metrics.csv.
%
%   See also FAST_TDMS_FOLDER2STRUCT, PLOT_RAPID_READOUT, RAPID_READOUT_PY,
%   RAPID_P1_PEAK, COMPUTE_M6_CENTROID, VECTOR_SUM_POLAR, FIND_PD_FROM_LUT.

    t_all = tic;
    here = fileparts(mfilename('fullpath'));
    addpath(here, fullfile(here, 'pipeline'), fullfile(here, '..', 'analyse_bar_DS'), fullfile(here, '..', 'helper'));

    session_folder = char(session_folder);
    if any(session_folder(end) == '/\'), session_folder = session_folder(1:end-1); end
    [~, sid] = fileparts(session_folder);
    opts = struct('Cell', sid, 'Out', fullfile(session_folder, 'rapid_readout'), 'Seed', [], 'Metrics', '', ...
                  'Figure', true, 'SavePng', true, 'Verbose', false, 'Stage1Test', false, 'Stage2Test', false, 'AgreeDeg', 10, 'FlagDeg', 5, ...
                  'Stage2', '', 'P1Peak', [], 'P1Summary', '');
    for k = 1:2:numel(varargin), opts.(varargin{k}) = varargin{k+1}; end
    if startsWith(opts.Out, '/Volumes/reiserlab')
        error('rapid_readout:ReadOnlyShare', 'the Janelia share is read-only: give an ''Out'' folder on the local disk');
    end
    if isempty(opts.Seed), opts.Seed = crc32(sid); end

    %% ---------------------------------------------------------------- read and analyze
    t0 = tic;
    [Log, info] = fast_tdms_folder2struct(session_folder, 'Verbose', opts.Verbose);
    t_read = toc(t0);
    t0 = tic;
    S = build_session(session_folder, sid, Log, session_folder);
    S.layout = info.layout;
    S2 = [];
    if ~isempty(opts.Stage2)
        st2 = char(opts.Stage2); if any(st2(end) == '/\'), st2 = st2(1:end-1); end
        t0b = tic; [Log2, info2] = fast_tdms_folder2struct(st2, 'Verbose', opts.Verbose); t_read = t_read + toc(t0b);
        cat_dir = session_folder;
        for c = {st2, fileparts(st2)}
            if isfolder(fullfile(c{1}, 'Patterns')) && isfolder(fullfile(c{1}, 'Functions')), cat_dir = c{1}; break; end
        end
        [~, sid2] = fileparts(st2);
        S2 = build_session(st2, sid2, Log2, cat_dir); S2.layout = info2.layout;
        assert(isequal(size(S2.fp), size(S.fp)) && max(abs(S2.fp(:) - S.fp(:))) == 0, 'rapid_readout:Stage2Pattern', 'the stage-2 FLASHES pattern differs from the stage-1 one');
        info.stage2 = info2; info.stage2_catalog_dir = cat_dir;
    end
    if opts.Stage1Test                                           % stage-1 emulation: no bar-flash trials
        keep = true(size(S.trials));
        for k = 1:numel(S.trials)
            if isKey(S.funcs, S.trials(k).fid), keep(k) = ~any(strcmp(S.funcs(S.trials(k).fid).cls, {'flash80', 'flash14'})); end
        end
        S.trials = S.trials(keep);
    end
    E = extract_events(S);
    if ~isempty(S2)
        E2 = extract_events(S2);
        if any(~cellfun(@isempty, E2.fl80)), E.fl80 = E2.fl80; end           % bar flashes from the stage-2 log
    end
    A = analyse(S, E);
    t_analysis = toc(t0);

    %% ---------------------------------------------------------------- placement PD, bar sets, geometry
    pp = placement_pd(A, opts.AgreeDeg, opts.FlagDeg);
    bs = bar_sets(A, S, pp.pd_placement);
    if opts.Stage2Test                                           % stage-2 emulation: only the 22 bars of the placement-frame sets were flashed
        keep = false(1, 88); keep([bs.perp_frames_1based, bs.par_frames_1based]) = true;
        A.flash80_amp995(~keep) = NaN; A.flash80_peak(~keep) = 0; A.flash80_under(~keep) = 0;
        A = tutl_axes(A, S); bs = bar_sets(A, S, pp.pd_placement);
    end
    n_bars = sum(isfinite(A.flash80_amp995)); A.tutl_available = isfinite(A.m6_pd) && isfinite(A.m6_or);
    rng(opts.Seed, 'twister');
    val = validate_metrics(A, opts, sid);
    tag = opts.Cell;
    if ~isfolder(opts.Out), mkdir(opts.Out); end
    fprintf('source: %s   (fast TDMS read %.2f s, analysis %.2f s; layout %s; seed %d)\n', session_folder, t_read, t_analysis, S.layout, opts.Seed);
    fprintf('cell %s  session %s  static Vm %.1f mV  snapped PD %g deg (reported metric pd_deg; margin of the 28 deg/s vector to the LUT boundary %.1f deg)\n', ...
        tag, sid, A.static_vm, A.pd, pp.margin_to_lut_boundary_deg);
    if ~bs.have_bars
        fprintf('read-out call 1: STAGE 1 session (sweeps only, no bar flashes) -> placement PD, bar sets, stage-2 bar list\n');
    else
        if isempty(S2), src2 = ''; else, src2 = sprintf(' [two-stage: sweeps from %s, bar flashes from %s]', session_folder, opts.Stage2); end
        fprintf('read-out call 2: bar flashes present (%d of 88 bars) -> center from the two placement-frame sets, window, sites, Prot3b and Protocol 4 tables%s\n', n_bars, src2);
    end
    if pp.pd_flag_gt5, flag = sprintf('  [pd_flag_gt5: the speeds disagree by more than %g deg]', opts.FlagDeg); else, flag = ''; end
    fprintf('placement PD %.1f deg: 28 deg/s vector %.1f deg (L = %.2f), 56 deg/s vector %.1f deg (L = %.2f), differ by %.1f deg -> rule %s%s\n', ...
        pp.pd_placement, pp.pd_vector_28, pp.vector_length_28, pp.pd_vector_56, pp.vector_length_56, pp.diff_28_56_deg, pp.pd_rule, flag);
    if ~pp.agree
        fprintf('   WARNING: the 28 and 56 deg/s vector PDs differ by %.1f deg (> %g deg); placement uses the 28 deg/s vector alone (28 deg/s %.1f deg, L %.2f; 56 deg/s %.1f deg, L %.2f)\n', ...
            pp.diff_28_56_deg, opts.AgreeDeg, pp.pd_vector_28, pp.vector_length_28, pp.pd_vector_56, pp.vector_length_56);
    end
    if bs.same_sets_as_tutl, same = ''; else, same = sprintf('  (TUTL snapped-PD choice: %d / %d)', A.g_pd, A.g_or); end
    fprintf('bar sets in the placement frame: perpendicular (PD-axis) set %d, parallel (orthogonal-axis) set %d%s\n', bs.g_perp, bs.g_par, same);

    if ~bs.have_bars
        %% ------------------------------------------------------------ stage 1: sweeps only -> stage-2 bar list and sequence
        [bars, seq2, total2] = stage2_tables(S, bs);
        fprintf('STAGE 1 session (no bar flashes). Stage 2: %d bars (frames of the 88-frame bar pattern) x 80 ms (+920) and 14 ms (+486), 3 repeats = %.1f min (5 repeats: %.1f min)\n', ...
            numel(bars), total2 / 60, total2 / 60 * 5 / 3);
        fprintf('   perpendicular set frames (along the PD axis): %s\n', mat2str(bs.perp_frames_1based));
        fprintf('   parallel set frames (along the orthogonal axis): %s\n', mat2str(bs.par_frames_1based));
        fprintf('sequencer inputs: usable_pd %d (%s)\n', pp.usable_pd, pp.usable_pd_rule);
        out = struct('cell', tag, 'session', sid, 'source', session_folder, 'stage', 1, 'seed', opts.Seed, 'pd_deg', A.pd, 'pd_vector_28', pp.pd_vector_28, 'pd_vector_56', pp.pd_vector_56, ...
            'pd_placement', pp.pd_placement, 'pd_rule', pp.pd_rule, 'pd_flag_gt5', pp.pd_flag_gt5, 'placement', pp, 'bar_sets', rmfield(bs, {'center_xy', 'e_pd', 'e_or'}), ...
            'sequencer_inputs', struct('usable_pd', pp.usable_pd, 'usable_pd_rule', pp.usable_pd_rule, 'pd_inputs', pp.pd_inputs, 'profiles_unimodal', 'n/a'), ...
            'stage2', struct('n_bars', numel(bars), 'reps', 3, 'total_s', round(total2, 1)), 'layout', S.layout, 'implementation', 'MATLAB rapid_readout.m');
        fid = fopen(fullfile(opts.Out, [tag '_stage1.json']), 'w'); fprintf(fid, '%s', jsonencode(out, 'PrettyPrint', true)); fclose(fid);
        write_csv(fullfile(opts.Out, [tag '_stage2_bars.csv']), bars, {'set', 'bar_orientation_deg', 'bar_pattern_frame_1based', 'rank_along_axis', 'col_center', 'row_center'}, {'bar_pattern_frame_1based', 'rank_along_axis'});
        write_csv(fullfile(opts.Out, [tag '_stage2_sequence.csv']), seq2, {'repeat', 'cls', 'order', 'bar_pattern_frame_1based', 'flash_ms', 'gray_ms', 't_start_s'}, {'repeat', 'order', 'bar_pattern_frame_1based', 'flash_ms', 'gray_ms'});
        fprintf('written: %s\n', fullfile(opts.Out, [tag '_{stage1.json,stage2_bars.csv,stage2_sequence.csv}']));
        R = struct('cell', tag, 'sid', sid, 'session_folder', session_folder, 'stage', 1, 'S', S, 'A', A, 'placement', pp, 'bar_sets', bs, 'stage2_bars', bars, 'stage2_sequence', seq2, ...
                   'validation', val, 'json', out, 'timing', struct('read_s', t_read, 'analysis_s', t_analysis, 'total_s', toc(t_all), 'read_info', info), 'out_dir', opts.Out);
        return;
    end

    %% ---------------------------------------------------------------- full session: geometry, sites, tables
    G = struct('C', bs.center_xy, 'e_pd', bs.e_pd, 'e_or', bs.e_or, 'screen', struct('col_min', 17, 'col_max', 192, 'row_min', 1, 'row_max', 48));
    G.center = to_arena(G, 0, 0);
    if A.tutl_available, C_tutl = tutl_center(A, S); else, C_tutl = [NaN NaN]; end
    [kept, grid] = lattice(G);
    [sites, notes] = place_sites(G);
    for i = 1:numel(sites)
        st = sites(i);
        if strcmp(st.region, 'outside')
            e = NaN;
        elseif strcmp(st.region, 'orthogonal')
            e = expected_response(bs, st.v_px, 'or');
        else
            e = expected_response(bs, st.s_px, 'pd');
        end
        if isfinite(e), sites(i).expected_ON_80ms_peak_mV = round(e, 1); else, sites(i).expected_ON_80ms_peak_mV = NaN; end
        q = find([kept.col_start] == st.col_start & [kept.row_start] == st.row_start, 1);
        sites(i).in_lattice = ~isempty(q);
        if isempty(q), sites(i).lattice_gen_frame_off = NaN; sites(i).lattice_gen_frame_on = NaN;
        else, sites(i).lattice_gen_frame_off = kept(q).gen_frame_off; sites(i).lattice_gen_frame_on = kept(q).gen_frame_on; end
    end
    [p3_frames, p3_blocks, p3_total] = prot3b_tables(G, sites);
    [p4_frames, p4_seq, p4_total, p4_counts, p4_item_s] = prot4_tables(G, kept);
    cc = G.center(1); cr = G.center(2);
    dist = struct('tutl_snapped_center', round(hypot(bs.center_xy(1) - C_tutl(1), bs.center_xy(2) - C_tutl(2)), 2));
    others = struct();
    p1 = opts.P1Peak;
    if isempty(p1) && ~isempty(opts.P1Summary) && isfile(opts.P1Summary)
        T1 = readtable(opts.P1Summary, 'TextType', 'char'); r1 = find(strcmp(T1.cell, tag), 1);
        if ~isempty(r1), p1 = [T1.ON_peak_x(r1), T1.ON_peak_y(r1)]; end
    end
    if ~isempty(p1), others.protocol1_ON_peak = p1(:)'; dist.protocol1_ON_peak = round(hypot(p1(1) - cc, p1(2) - cr), 2); end
    pmp = profile_metrics(bs.perp_pos_px, bs.perp_peak); pmo = profile_metrics(bs.par_pos_px, bs.par_peak);
    uni = struct('pd', halfmax_peaks(bs.perp_peak), 'or', halfmax_peaks(bs.par_peak));
    seq_in = struct('usable_pd', pp.usable_pd, 'usable_pd_rule', pp.usable_pd_rule, 'pd_inputs', pp.pd_inputs, ...
                    'profiles_unimodal', 'n/a', 'profiles_unimodal_note', 'no justified rule: counting local maxima at or above half of the profile maximum calls 18 of 25 good TUTL cells multimodal (broad orthogonal plateaus); the inputs are given instead', ...
                    'profile_halfmax_peaks', uni, 'profile_pd_pos_px', bs.perp_pos_px, 'profile_pd_peak_mV', round(bs.perp_peak, 2), 'profile_or_pos_px', bs.par_pos_px, 'profile_or_peak_mV', round(bs.par_peak, 2), ...
                    'centroid_offset_placement_px', struct('s', round(pmp.centroid, 2), 'v', round(pmo.centroid, 2)), ...
                    'width68_placement_px', struct('s', round(pmp.width68, 2), 'v', round(pmo.width68, 2)), ...
                    'center_expected_ON_80ms_peak_mV', sites(1).expected_ON_80ms_peak_mV);

    %% ---------------------------------------------------------------- print
    if A.tutl_available, tutl_txt = sprintf('TUTL sets: %d / %d', A.m6_pd, A.m6_or); else, tutl_txt = 'TUTL sets not presented in this session'; end
    fprintf('RF center (placement frame, M6 x M6): col %.1f, row %.1f (1-based); M6 bars %d / %d of 11 (%s)   centroid offset s %+.1f, v %+.1f px\n', ...
        cc, cr, bs.perp_m6_1based, bs.par_m6_1based, tutl_txt, A.centroid_pd, A.centroid_or);
    if A.tutl_available, fprintf('   distance to tutl_snapped_center: %.1f px\n', dist.tutl_snapped_center); end
    if isfield(dist, 'protocol1_ON_peak'), fprintf('   distance to protocol1_ON_peak (col %.1f row %.1f): %.1f px\n', p1(1), p1(2), dist.protocol1_ON_peak); end
    fprintf('sequencer inputs: usable_pd %d (%s); half-max peaks of the 80 ms profiles PD %d, orth %d (no unimodality verdict); placement-frame centroid offset s %+.1f, v %+.1f px, width68 s %.1f, v %.1f px; center expected ON 80 ms peak %s mV\n', ...
        seq_in.usable_pd, pp.usable_pd_rule, uni.pd, uni.or, pmp.centroid, pmo.centroid, pmp.width68, pmo.width68, num2str(seq_in.center_expected_ON_80ms_peak_mV));
    if isempty(val)
        fprintf('validation vs metrics.csv (TUTL quantities): not requested\n');
    else
        fprintf('validation vs metrics.csv (TUTL quantities): %s %s\n', val.status, val.detail);
    end
    fprintf('PROT3B  26 blocks, %.1f min; pattern %d frames (background, %d sites x OFF/ON at 4 px, %d 3-px OFF squares); AO levels -1/-0.5/0/+0.5/+1 pA = %s V\n', ...
        p3_total / 60, numel(p3_frames), numel(sites), sum(strcmp({p3_frames.item}, 'mi9_off')), strjoin(arrayfun(@(c) sprintf('%+.3f', c * 0.25), [-1 -.5 0 .5 1], 'UniformOutput', false), '/'));
    for k = 1:numel(notes), fprintf('   WARNING: %s\n', notes{k}); end
    fprintf('   sites (4-px squares; col/row = 1-based start pixel; frames = Prot3b pattern frame numbers OFF/ON):\n');
    for i = 1:numel(sites)
        st = sites(i); fo = 1 + i;
        if isfinite(st.expected_ON_80ms_peak_mV), ex = sprintf('%.1f mV', st.expected_ON_80ms_peak_mV); else, ex = 'n/a'; end
        if isempty(st.note), nt = ''; else, nt = ['  [' st.note ']']; end
        fprintf('   %-10s target s %+3d v %+3d -> s %+5.1f v %+5.1f px  start col %3d row %3d  frames %2d/%2d  expected ON 80 ms peak %s%s\n', ...
            st.region, st.s_target, st.v_target, st.s_px, st.v_px, st.col_start, st.row_start, fo, fo + numel(sites), ex, nt);
    end
    mi9 = [p3_frames(strcmp({p3_frames.item}, 'mi9_off')).frame];
    fprintf('   Mi9 holding series: 3-px OFF squares at s = [0, -4, -8] (frames %s), 500 ms, 1-s period, 6 repeats per block, blocks 21-26 at [-1, 0, 1, 1, 0, -1] pA\n', mat2str(mi9));
    fprintf('PROTOCOL 4  %.1f min, 6 passes (lattice in passes 1, 3, 5); pattern %d frames: lattice %d (%d positions x OFF/ON, 4-px at 2-px steps, px_rng rows %d-%d cols %d-%d, %d off-screen dropped), core %d (3-px ON, 1-px steps), line %d (2-px ON), off3 %d (3-px OFF); per item lattice %.1f min, core %.1f min, line %.1f min, off3 %.1f min\n', ...
        p4_total / 60, numel(p4_frames), p4_counts.lattice, grid.n_kept, grid.px_rng.row_start, grid.px_rng.row_end, grid.px_rng.col_start, grid.px_rng.col_end, ...
        grid.n_dropped_off_screen, p4_counts.core, p4_counts.line, p4_counts.off3, p4_item_s.lattice / 60, p4_item_s.core / 60, p4_item_s.line / 60, p4_item_s.off3 / 60);

    %% ---------------------------------------------------------------- files
    out = struct('cell', tag, 'session', sid, 'source', session_folder, 'stage', 2, 'analysis_seconds', round(t_read + t_analysis, 1), 'seed', opts.Seed, ...
        'pd_deg', A.pd, 'pd_vector_28', pp.pd_vector_28, 'pd_vector_56', pp.pd_vector_56, 'pd_placement', pp.pd_placement, 'pd_rule', pp.pd_rule, 'pd_flag_gt5', pp.pd_flag_gt5, 'placement', pp, ...
        'bar_sets', rmfield(bs, {'center_xy', 'e_pd', 'e_or', 'perp_peak', 'par_peak', 'perp_under', 'par_under'}), 'm6_pd_1based_tutl', A.m6_pd, 'm6_or_1based_tutl', A.m6_or, ...
        'static_vm_mV', round(double(A.static_vm), 2), 'rf_center_col_row_1based', [round(cc, 2), round(cr, 2)], ...
        'centroid_offset_px', struct('s', round(A.centroid_pd, 2), 'v', round(A.centroid_or, 2)), ...
        'window_px', struct('s_lead', -12, 's_trail', 16, 'v_half', 8), 'core_px', struct('s_half', 6, 'v_half', 4), 'lattice', grid, 'deg_per_px', 1.25, ...
        'sites', sites, 'site_warnings', {notes}, ...
        'prot3b', struct('n_blocks', numel(p3_blocks), 'total_s', round(p3_total, 1), 'timing', struct('flash_ms', 160, 'gray_ms', 440, 'settle_s', 3.0, 'hold_flash_ms', 500, 'hold_period_ms', 1000, 'hold_reps', 6), ...
                         'order_blocks_1_20_pA', p3b_order(), 'order_blocks_21_26_pA', [-1 0 1 1 0 -1], 'ao_volts_per_pA', 0.25, 'n_frames', numel(p3_frames)), ...
        'protocol4', struct('total_s', round(p4_total, 1), 'counts', p4_counts, 'seconds_per_item', p4_item_s, 'n_frames', numel(p4_frames)), ...
        'other_centers', others, 'distance_px', dist, 'sequencer_inputs', seq_in, 'validation', val, 'layout', S.layout, 'implementation', 'MATLAB rapid_readout.m', ...
        'conventions', 's + = toward placement PD = trailing/proximal; v + = 90 deg clockwise from PD; arena col/row 1-based, row 1 = top; frame 1 = background in both patterns');
    fid = fopen(fullfile(opts.Out, [tag '_locations.json']), 'w'); fprintf(fid, '%s', jsonencode(out, 'PrettyPrint', true)); fclose(fid);
    write_csv(fullfile(opts.Out, [tag '_lattice.csv']), kept, {'size_px', 'col_start', 'row_start', 'col_center', 'row_center', 's_px', 'v_px', 's_target', 'v_target', 'on_screen', 'gen_frame_off', 'gen_frame_on'}, ...
        {'size_px', 'col_start', 'row_start', 'gen_frame_off', 'gen_frame_on'});
    write_csv(fullfile(opts.Out, [tag '_sites.csv']), sites, {'size_px', 'col_start', 'row_start', 'col_center', 'row_center', 's_px', 'v_px', 's_target', 'v_target', 'on_screen', 'region', 'note', 'expected_ON_80ms_peak_mV', 'in_lattice', 'lattice_gen_frame_off', 'lattice_gen_frame_on'}, ...
        {'size_px', 'col_start', 'row_start', 's_target', 'v_target', 'lattice_gen_frame_off', 'lattice_gen_frame_on'});
    write_csv(fullfile(opts.Out, [tag '_prot3b_frames.csv']), p3_frames, {'frame', 'item', 'size_px', 'polarity', 'col_start', 'row_start', 'col_center', 'row_center', 's_px', 'v_px', 's_target', 'v_target', 'region'}, ...
        {'frame', 'size_px', 'col_start', 'row_start', 's_target', 'v_target'});
    write_csv(fullfile(opts.Out, [tag '_prot3b_blocks.csv']), p3_blocks, {'block', 'phase', 'current_pA', 'ao_volts', 'settle_s', 'n_flashes', 'flash_ms', 'gray_ms', 'duration_s', 't_start_s', 'frames'}, ...
        {'block', 'n_flashes', 'flash_ms', 'gray_ms'});
    write_csv(fullfile(opts.Out, [tag '_prot4_frames.csv']), p4_frames, {'frame', 'item', 'size_px', 'polarity', 'col_start', 'row_start', 'col_center', 'row_center', 's_px', 'v_px'}, ...
        {'frame', 'size_px', 'col_start', 'row_start'});
    write_csv(fullfile(opts.Out, [tag '_prot4_sequence.csv']), p4_seq, {'pass', 'order', 'frame', 'item', 'size_px', 'polarity', 'col_start', 'row_start', 's_px', 'v_px', 'flash_ms', 'gray_ms', 't_start_s'}, ...
        {'pass', 'order', 'frame', 'size_px', 'col_start', 'row_start', 'flash_ms', 'gray_ms'});

    %% ---------------------------------------------------------------- figure
    fig = [];
    if opts.Figure
        P = struct('cell', tag, 'pd', pp.pd_placement, 'pd_vector', pp.pd_vector_28, 'pd_vector_56', pp.pd_vector_56, 'pd_snapped', A.pd, 'center', G.center, 'center_xy', G.C, 'e_pd', G.e_pd, 'e_or', G.e_or, ...
                   'bar_cx', S.cx + 1, 'bar_cy', S.cy + 1, 'g_pd', bs.g_perp, 'g_or', bs.g_par, 'lattice', kept, 'sites', sites, 'p4_frames', p4_frames, ...
                   'window', struct('s_lead', -12, 's_trail', 16, 'v_half', 8));
        if ~isempty(p1), P.others = struct('name', 'protocol1_ON_peak', 'col', p1(1), 'row', p1(2)); end
        P.profiles.pd = struct('pos', bs.perp_pos_px, 'peak', bs.perp_peak, 'under', bs.perp_under, 'm6', bs.perp_m6_1based, 'set', bs.g_perp);
        P.profiles.or = struct('pos', bs.par_pos_px, 'peak', bs.par_peak, 'under', bs.par_under, 'm6', bs.par_m6_1based, 'set', bs.g_par);
        fig = plot_rapid_readout(P);
        if opts.SavePng
            exportgraphics(fig, fullfile(opts.Out, [tag '_locations.png']), 'Resolution', 110);
        end
    end
    t_total = toc(t_all);
    fprintf('written: %s   (total %.1f s: read %.2f, analysis %.2f, tables+files+figure %.2f)\n', fullfile(opts.Out, [tag '_{locations.json,sites.csv,lattice.csv,prot3b_frames.csv,prot3b_blocks.csv,prot4_frames.csv,prot4_sequence.csv,locations.png}']), ...
        t_total, t_read, t_analysis, t_total - t_read - t_analysis);

    R = struct('cell', tag, 'sid', sid, 'session_folder', session_folder, 'stage', 2, 'S', S, 'A', A, 'placement', pp, 'bar_sets', bs, 'G', G, 'sites', sites, 'site_warnings', {notes}, 'lattice', kept, 'grid', grid, ...
               'prot3b_frames', p3_frames, 'prot3b_blocks', p3_blocks, 'prot4_frames', p4_frames, 'prot4_sequence', p4_seq, 'validation', val, 'json', out, ...
               'timing', struct('read_s', t_read, 'analysis_s', t_analysis, 'total_s', t_total, 'read_info', info), 'figure', fig, 'out_dir', opts.Out);
end


%% ===================================================================== session structure
function S = build_session(sdir, sid, Log, cat_dir)
% Mirror of pharm_protocol.load_session: the frame channel, Vm (mV, single), the trial table from the command log, the
% pattern/function catalogs from the file names and the bar geometry from the FLASHES pattern.
    S = struct('sid', sid, 'sdir', sdir);
    % catalogs
    if nargin < 4, cat_dir = sdir; end
    pf = dir(fullfile(cat_dir, 'Patterns', '*.mat')); pats = containers.Map('KeyType', 'double', 'ValueType', 'any');
    for k = 1:numel(pf)
        n = pf(k).name; id = str2double(n(1:4));
        t = regexp(n, '_(\d)-8pi_', 'tokens', 'once');
        if ~isempty(t), orient = mod(90 - 22.5 * str2double(t{1}), 360);
        elseif ~isempty(regexp(n, '_pi_', 'once')), orient = 90;
        else, orient = NaN; end
        pats(id) = struct('name', n, 'path', fullfile(pf(k).folder, n), 'orient', orient);
    end
    ff = dir(fullfile(cat_dir, 'Functions', '*.mat')); funcs = containers.Map('KeyType', 'double', 'ValueType', 'any');
    for k = 1:numel(ff)
        n = ff(k).name; id = str2double(n(1:4));
        t = regexp(n, '(\d+)dps', 'tokens', 'once');
        if ~isempty(t), cls = 'sweep'; speed = str2double(t{1});
        elseif contains(n, 'flashes_80ms'), cls = 'flash80'; speed = NaN;
        elseif contains(n, 'flashes_14ms'), cls = 'flash14'; speed = NaN;
        elseif contains(n, 'static'), cls = 'static'; speed = NaN;
        else, cls = 'other'; speed = NaN; end
        funcs(id) = struct('name', n, 'cls', cls, 'speed', speed, 'flip', contains(upper(n), 'FLIP'));
    end
    S.pats = pats; S.funcs = funcs;
    % channels
    av = Log.ADC.Volts; at = double(Log.ADC.Time);
    v0 = av(1, :);
    if mean(abs(v0 - round(v0)) < 1e-3) > 0.9
        fr = int16(round(v0)); assert(max(abs(v0 - double(fr))) < 1e-3);
        vm = single(av(2, :) * 10); ticks = at(2, :); S.i_monitor_pa = NaN; S.layout = 'ADC0 frames / ADC1 Vm';
    else
        pos = double(Log.Frames.Position(1, :)); ft = double(Log.Frames.Time(1, :));
        ticks = at(1, :);
        idx = discretize(ticks, [ft, Inf]);                      % number of frame-log entries with time <= tick (0 before the first)
        fr = zeros(1, numel(ticks), 'int16'); ok = ~isnan(idx); fr(ok) = int16(pos(idx(ok)));
        vm = single(av(1, :) * 10); S.i_monitor_pa = median(av(2, :)) * 2000; S.layout = 'ADC1 Vm / ADC2 current + frame log';
    end
    S.fr = fr; S.vm = vm; S.n = numel(fr);
    % trials from the command log (0-based sample index of Start-Display; end = next start)
    names = Log.Commands.Name; data = Log.Commands.Data; ct = double(Log.Commands.Time);
    pid = NaN; fidn = NaN; sd_t = []; sd_pid = []; sd_fid = [];
    for k = 1:numel(names)
        switch names{k}
            case 'Set Pattern ID', pid = hex2dec(strtok(data{k}, '/'));
            case 'Set Pattern Function ID', fidn = hex2dec(strtok(data{k}, '/'));
            case 'Start-Display'
                sd_t(end+1) = ct(k); sd_pid(end+1) = pid; sd_fid(end+1) = fidn; %#ok<AGROW>
        end
    end
    c = discretize(sd_t, [ticks, Inf]);                            % count of ticks <= t (NaN before the first tick)
    c(isnan(c)) = 0;
    hit = c >= 1; hit(hit) = ticks(c(hit)) == sd_t(hit); c(hit) = c(hit) - 1;   % searchsorted 'left': count of ticks < t
    trials = struct('start', num2cell(c), 'end', num2cell([c(2:end), S.n]), 'pid', num2cell(sd_pid), 'fid', num2cell(sd_fid));
    S.trials = trials;
    % FLASHES pattern geometry
    keys = cell2mat(pats.keys); fl = keys(cellfun(@(k) contains(pats(k).name, 'FLASHES'), num2cell(keys)));
    assert(numel(fl) == 1, 'rapid_readout:Flashes', 'expected one FLASHES pattern, found %d', numel(fl));
    pz = load(pats(fl).path, 'pattern'); Pm = single(pz.pattern.Pats);
    is_on = max(Pm(:)) > 4; S.is_on = is_on; S.flash_pattern = pats(fl).name;
    if is_on, W = min(max((Pm - 4) / 11, 0), 1); else, W = min(max((4 - Pm) / 4, 0), 1); end
    [H, Wd, ~] = size(W);
    fp = double(reshape(W(:, :, 2:89), H * Wd, 88))';            % 88 x (H*W), column-major pixel order
    [rr, cc] = ndgrid(0:H-1, 0:Wd-1); X = cc(:)'; Y = rr(:)';
    a = sum(fp, 2); S.cx = (fp * X')' ./ a'; S.cy = (fp * Y')' ./ a';
    go = zeros(1, 8);
    for g = 0:7
        pts = [];
        for k = 1:11
            m = fp(g * 11 + k, :) > 0.5; pts = [pts; X(m)', -Y(m)']; %#ok<AGROW>
        end
        pts = pts - mean(pts, 1);
        [~, ~, V] = svd(pts, 'econ'); u = V(:, 1);
        go(g + 1) = mod(atan2d(u(2), u(1)), 180);
    end
    S.g_orient = go; S.fp = fp; S.shape = [H, Wd];
end


%% ===================================================================== events (mirror of pharm_protocol.events, PD and center parts)
function E = extract_events(S)
    fr = S.fr; vm = S.vm; n = S.n;
    E.static = []; E.sweep_keys = zeros(0, 2);
    E.sw28 = containers.Map('KeyType', 'char', 'ValueType', 'any');     % direction -> list of 28 dps trials
    E.sw56 = containers.Map('KeyType', 'char', 'ValueType', 'any');     % direction -> list of 56 dps trials (1 kHz traces, for the placement PD)
    E.fl80 = repmat({{}}, 1, 88);                                        % bar -> list of flashes
    for t = S.trials
        if ~isKey(S.funcs, t.fid), continue; end
        f = S.funcs(t.fid);
        if isKey(S.pats, t.pid), p = S.pats(t.pid); else, p = []; end
        seg = fr(t.start + 1:t.end);
        switch f.cls
            case 'static'
                E.static(end+1) = median(vm(t.start + 1:t.end));
            case 'sweep'
                zero = find(seg == 0);
                if isempty(zero), continue; end
                after = find(seg(zero(1):end) ~= 0);
                if isempty(after), continue; end
                on = t.start + zero(1) - 1 + after(1) - 1;              % 0-based motion onset
                nz = find(seg ~= 0); off = t.start + nz(end);           % 0-based exclusive motion end
                if off - on < 500, continue; end
                if isempty(p) || isnan(p.orient), continue; end
                direction = mod(p.orient - 90 + 180 * f.flip, 360);
                a = on - 3000; b = off + 3000;
                if a < 0 || b > n, continue; end
                E.sweep_keys(end+1, :) = [f.speed, direction];
                if f.speed == 56
                    snip = double(vm(a + 1:b)); base = snip(1:3000);
                    key = sprintf('%.4f', direction);
                    item = struct('tr', decimate10(snip - mean(base)), 'motion_ms', (off - on) / 10);
                    if isKey(E.sw56, key), lst = E.sw56(key); lst{end+1} = item; E.sw56(key) = lst; else, E.sw56(key) = {item}; end
                end
                if f.speed == 28
                    r0 = max(on - 10000, t.start);
                    key = sprintf('%.4f', direction);
                    item = struct('raw', vm(r0 + 1:t.end), 'raw_on', on - r0, 'motion_n', off - on);
                    if isKey(E.sw28, key), lst = E.sw28(key); lst{end+1} = item; E.sw28(key) = lst; else, E.sw28(key) = {item}; end
                end
            case 'flash80'
                ch = find(diff(seg) ~= 0);                               % 0-based positions of run starts (Python: flatnonzero(diff)+1)
                bnd = [0, ch, numel(seg)];
                for j = 1:numel(bnd) - 1
                    x = bnd(j); y = bnd(j + 1); v = double(seg(x + 1));
                    if v == 0 || x == 0 || (y - x) > 1600, continue; end
                    on = t.start + x; a = on - 2500; b = on + 9000;
                    if a < 0 || b > n, continue; end
                    snip = double(vm(a + 1:b)); base = snip(1:2500);
                    tr = decimate10(snip - mean(base));
                    if on >= 5000 && on + 1551 <= n, raw = vm(on - 5000 + 1:on + 1551); else, raw = []; end
                    E.fl80{v}{end+1} = struct('tr', tr, 'raw', raw);
                end
        end
    end
end

function y = decimate10(x)
    m = floor(numel(x) / 10) * 10; y = mean(reshape(x(1:m), 10, m / 10), 1);
end

function y = smooth5(x)
    y = conv(x, ones(1, 5) / 5, 'same');
end


%% ===================================================================== analysis (mirror of pharm_protocol.analyse_epoch / analyse_group)
function A = analyse(S, E)
    if isempty(E.static), A.static_vm = NaN; else, A.static_vm = median(E.static); end
    % direction list: from the 56 dps sweeps as in analyse_group, else from the 28 dps sweeps
    k56 = E.sweep_keys(E.sweep_keys(:, 1) == 56, 2);
    if isempty(k56), k56 = E.sweep_keys(E.sweep_keys(:, 1) == 28, 2); end
    dirs = sort(unique(k56))';
    % PD amplitude per direction (compute_bar_sweep_responses convention on the repeat mean at 10 kHz)
    amp28 = nan(size(dirs));
    for i = 1:numel(dirs)
        key = sprintf('%.4f', dirs(i));
        if ~isKey(E.sw28, key), continue; end
        lst = E.sw28(key);
        pre_n = min(cellfun(@(x) x.raw_on, lst)); post_n = min(cellfun(@(x) numel(x.raw) - x.raw_on, lst));
        M = zeros(numel(lst), pre_n + post_n);
        for r = 1:numel(lst), M(r, :) = double(lst{r}.raw(lst{r}.raw_on - pre_n + 1:lst{r}.raw_on + post_n)); end
        m = mean(M, 1); o = pre_n;
        if o > 1000, bl = mean(m(max(o - 9000, 0) + 1:o - 1000)); else, bl = mean(m(1:o)); end
        mot_n = min(cellfun(@(x) x.motion_n, lst));
        amp28(i) = abs(np_percentile(m(o + 1:min(o + mot_n + 1500, numel(m))), 98) - bl);
    end
    amp56 = nan(size(dirs));
    for i = 1:numel(dirs)
        key = sprintf('%.4f', dirs(i));
        if ~isKey(E.sw56, key), continue; end
        lst = E.sw56(key); L = min(cellfun(@(x) numel(x.tr), lst)); M = zeros(numel(lst), L);
        for r = 1:numel(lst), M(r, :) = lst{r}.tr(1:L); end
        m = mean(M, 1); mot = floor(median(cellfun(@(x) x.motion_ms, lst))); sm = smooth5(m);
        w0 = 300; w1 = min(w0 + mot + 150, numel(m));                 % [motion onset, motion end + 150 ms) at 1 kHz
        amp56(i) = max(sm(w0 + 1:w1));
    end
    A.dirs = dirs; A.amp28 = amp28; A.amp56 = amp56;
    ok = isfinite(amp28);
    [~, ang] = vector_sum_polar(amp28(ok), deg2rad(dirs(ok)));
    ang = mod(rad2deg(ang), 360);
    [~, i] = min(abs(wrap180(dirs - ang)));
    A.pd = dirs(i); A.pd_vector = ang; A.nd = mod(A.pd + 180, 360);
    % 80 ms bar-flash amplitudes (compute_pos_peak_amplitudes convention) and the 1 kHz profiles
    a995 = nan(1, 88); tr = nan(88, 1150);
    for v = 1:88
        lst = E.fl80{v}; if isempty(lst), continue; end
        raws = lst(cellfun(@(x) ~isempty(x.raw), lst));
        if ~isempty(raws)
            M = zeros(numel(raws), 6551);
            for r = 1:numel(raws), M(r, :) = double(raws{r}.raw); end
            m = mean(M, 1); a995(v) = np_percentile(m(5001:6551) - mean(m(1:5000)), 99.5);
        end
        L = min(cellfun(@(x) numel(x.tr), lst)); M = zeros(numel(lst), L);
        for r = 1:numel(lst), M(r, :) = lst{r}.tr(1:L); end
        tr(v, 1:L) = mean(M, 1);
    end
    tr(isnan(tr)) = 0; sm = zeros(size(tr));
    for v = 1:88, sm(v, :) = smooth5(tr(v, :)); end
    peak = max(sm(:, 251:550), [], 2)'; under = min(sm(:, 401:750), [], 2)';
    A.flash80_amp995 = a995; A.flash80_peak = peak; A.flash80_under = under;
    A = tutl_axes(A, S);
end

function A = tutl_axes(A, S)
% The TUTL quantities: bar sets chosen from the snapped PD and their M6 centers (NaN when a set has no bar flashes, e.g. a stage-2 session).
    go = S.g_orient; a995 = A.flash80_amp995; peak = A.flash80_peak; under = A.flash80_under;
    [~, A.g_pd] = min(abs(mod(go - (A.pd + 90) + 90, 180) - 90)); A.g_pd = A.g_pd - 1;     % bars perpendicular to PD -> PD-axis profile
    [~, A.g_or] = min(abs(mod(go - A.pd + 90, 180) - 90)); A.g_or = A.g_or - 1;            % bars parallel to PD -> orthogonal-axis profile
    e_pd = [cosd(A.pd), sind(A.pd)]; e_or = [cosd(A.pd - 90), sind(A.pd - 90)];
    A.e_pd = e_pd; A.e_or = e_or;
    axn = {'pd', 'or'}; gs = [A.g_pd, A.g_or]; es = {e_pd, e_or};
    for j = 1:2
        g = gs(j); idx = g * 11 + (1:11);
        pos = [S.cx(idx)', -S.cy(idx)'] * es{j}';
        [~, order] = sort(pos);                                          % ascending along the axis = the paper's pos_order
        amp = a995(idx(order));
        if all(isnan(amp))
            A.(['bar_pos_' axn{j}]) = (pos - pos(order(6)))'; A.(['kc_' axn{j}]) = NaN; A.(['m6_' axn{j}]) = NaN;
            A.(['profile_' axn{j}]) = nan(1, 11); A.(['under_' axn{j}]) = nan(1, 11); A.(['amp995_' axn{j}]) = amp(order);
            A.(['centroid_' axn{j}]) = NaN; A.(['width68_' axn{j}]) = NaN;
            continue;
        end
        amp(isnan(amp)) = 0;
        m6 = compute_m6_centroid(max(amp, 0), 0.68);
        kc_sorted = m6.centroid_int; kc = order(kc_sorted);              % center bar, 1-based within the set
        A.(['bar_pos_' axn{j}]) = (pos - pos(kc))'; A.(['kc_' axn{j}]) = kc - 1; A.(['m6_' axn{j}]) = kc_sorted;
        A.(['profile_' axn{j}]) = peak(idx); A.(['under_' axn{j}]) = under(idx); A.(['amp995_' axn{j}]) = a995(idx);
        pm = profile_metrics(A.(['bar_pos_' axn{j}]), peak(idx));
        A.(['centroid_' axn{j}]) = pm.centroid; A.(['width68_' axn{j}]) = pm.width68;
    end
end

function pm = profile_metrics(pos, prof)
    w = prof; w(isnan(w)) = 0; w = max(w, 0);
    if sum(w) <= 0, pm = struct('centroid', NaN, 'width68', NaN); return; end
    pm.centroid = sum(w .* pos) / sum(w);
    [~, k] = max(w); lo = k; hi = k;
    while sum(w(lo:hi)) < 0.68 * sum(w)
        if lo > 1, left = w(lo - 1); else, left = -1; end
        if hi < numel(w), right = w(hi + 1); else, right = -1; end
        if left >= right, lo = lo - 1; else, hi = hi + 1; end
    end
    pm.width68 = abs(pos(hi) - pos(lo)) + abs(pos(2) - pos(1));
end

function y = np_percentile(x, q)
% numpy.percentile default ('linear'): index = q/100 * (n-1), linear interpolation between order statistics.
    x = sort(x(:)); n = numel(x); r = q / 100 * (n - 1); i = floor(r); f = r - i;
    if i + 1 >= n, y = x(n); else, y = x(i + 1) + f * (x(i + 2) - x(i + 1)); end
end

function d = wrap180(d)
    d = mod(d + 180, 360) - 180;
end

function e = expected_response(bs, x, axn)
% 80 ms bar peak (mV) interpolated along the placement PD axis (perpendicular set) at s, or the orthogonal axis (parallel set) at v.
    if strcmp(axn, 'pd'), pos = bs.perp_pos_px; prof = bs.perp_peak; else, pos = bs.par_pos_px; prof = bs.par_peak; end
    [p, o] = sort(pos); e = interp1(p, prof(o), x, 'linear', NaN);
end


%% ===================================================================== placement PD and bar sets (mirror of rapid_prot2_locations.py placement_pd / bar_sets)
function [ang, len] = vec_pd(dirs, amp)
    ok = isfinite(amp) & amp > 0;
    v = sum(amp(ok) .* exp(1i * deg2rad(dirs(ok)))); ang = mod(rad2deg(angle(v)), 360); len = abs(v) / sum(amp(ok));
end

function pp = placement_pd(A, agree_deg, flag_deg)
    [v28, l28] = vec_pd(A.dirs, A.amp28); [v56, l56] = vec_pd(A.dirs, A.amp56);
    d = abs(wrap180(v28 - v56));
    if d <= agree_deg
        place = mod(rad2deg(angle(exp(1i * deg2rad(v28)) + exp(1i * deg2rad(v56)))), 360); rule = 'mean';
    else
        place = v28; rule = '28_only';
    end
    assert(abs(wrap180(v28 - A.pd_vector)) < 0.051, 'rapid_readout:Placement', 'the 28 deg/s vector differs from analyse_group''s pd_vector');
    pp = struct('pd_deg', A.pd, 'pd_vector_28', round(v28, 2), 'pd_vector_56', round(v56, 2), 'vector_length_28', round(l28, 3), 'vector_length_56', round(l56, 3), ...
                'diff_28_56_deg', round(d, 2), 'pd_placement', round(place, 2), 'pd_rule', rule, 'pd_flag_gt5', d > flag_deg, 'agree', d <= agree_deg, 'agree_deg', agree_deg, 'flag_deg', flag_deg, ...
                'snapped_minus_vector28_deg', round(abs(wrap180(A.pd - v28)), 2), 'snapped_minus_placement_deg', round(abs(wrap180(A.pd - place)), 2), ...
                'margin_to_lut_boundary_deg', round(11.25 - abs(wrap180(v28 - A.pd)), 2));
    % usable_pd: the 28 deg/s tuning can carry a vector PD when every direction has an amplitude and the normalized vector length is
    % at least 0.1 (the sequencer's own floor; the 25 TUTL sessions span 0.14-0.48)
    n_ok = sum(isfinite(A.amp28)); [~, ipd] = min(abs(wrap180(A.dirs - A.pd))); [~, ind] = min(abs(wrap180(A.dirs - A.nd)));
    pp.usable_pd = n_ok == 16 && l28 >= 0.1;
    pp.usable_pd_rule = 'all 16 directions measured and normalized 28 deg/s vector length >= 0.1';
    pp.pd_inputs = struct('vector_length_28', round(l28, 3), 'vector_length_56', round(l56, 3), 'n_directions_28', n_ok, ...
                          'amp28_max_mV', round(max(A.amp28), 2), 'amp28_at_pd_mV', round(A.amp28(ipd), 2), 'amp28_at_nd_mV', round(A.amp28(ind), 2));
end

function n = halfmax_peaks(prof)
% number of local maxima of the positive part of an 11-point profile at or above half of its maximum
    w = prof; w(isnan(w)) = 0; w = max(w, 0); n = 0;
    if max(w) <= 0, return; end
    for k = 1:numel(w)
        left = k == 1 || w(k) > w(k - 1); right = k == numel(w) || w(k) >= w(k + 1);
        if left && right && w(k) >= 0.5 * max(w), n = n + 1; end
    end
end

function bs = bar_sets(A, S, pd_place)
% The two bar sets of the placement frame: perpendicular bars (orientation nearest placement + 90) give the PD-axis profile, parallel
% bars (nearest the placement PD) the orthogonal profile; per set the bars ordered along the axis and, with bar flashes, the M6 center.
    go = S.g_orient; e_pd = [cosd(pd_place), sind(pd_place)]; e_or = [cosd(pd_place - 90), sind(pd_place - 90)];
    [~, g_perp] = min(abs(mod(go - (pd_place + 90) + 90, 180) - 90)); g_perp = g_perp - 1;
    [~, g_par] = min(abs(mod(go - pd_place + 90, 180) - 90)); g_par = g_par - 1;
    amp = A.flash80_amp995; have = any(isfinite(amp));
    bs = struct('g_perp', g_perp, 'g_par', g_par, 'tutl_g_pd', A.g_pd, 'tutl_g_or', A.g_or, 'same_sets_as_tutl', g_perp == A.g_pd && g_par == A.g_or, 'have_bars', have);
    names = {'perp', 'par'}; gs = [g_perp, g_par]; es = {e_pd, e_or}; kcs = [NaN NaN];
    for j = 1:2
        g = gs(j); idx = g * 11 + (1:11);
        pos = ([S.cx(idx)', -S.cy(idx)'] * es{j}')'; [~, order] = sort(pos);
        if have
            a = amp(idx(order)); a(isnan(a)) = 0; m6 = compute_m6_centroid(max(a, 0), 0.68);
            kc = order(m6.centroid_int); kcs(j) = kc;
            bs.([names{j} '_m6_1based']) = m6.centroid_int; bs.([names{j} '_kc']) = kc - 1; bs.([names{j} '_pos_px']) = round(pos - pos(kc), 2);
            bs.([names{j} '_peak']) = A.flash80_peak(idx); bs.([names{j} '_under']) = A.flash80_under(idx);
        end
        bs.([names{j} '_frames_1based']) = idx(order); bs.([names{j} '_order']) = order - 1;
    end
    bs.e_pd = e_pd; bs.e_or = e_or; bs.center_xy = [];
    if have
        kp = g_perp * 11 + kcs(1); ko = g_par * 11 + kcs(2);
        r_pd = [S.cx(kp), -S.cy(kp)]; r_or = [S.cx(ko), -S.cy(ko)];
        bs.center_xy = (r_pd * e_pd') * e_pd + (r_or * e_or') * e_or;
    end
end

function C = tutl_center(A, S)
% 2-D center from the TUTL sets and the snapped PD (the rule before 26 Oct 15:05), kept for the distance check.
    k_pd = A.g_pd * 11 + A.kc_pd + 1; k_or = A.g_or * 11 + A.kc_or + 1;
    r_pd = [S.cx(k_pd), -S.cy(k_pd)]; r_or = [S.cx(k_or), -S.cy(k_or)];
    C = (r_pd * A.e_pd') * A.e_pd + (r_or * A.e_or') * A.e_or;
end

function [rows, seq, t] = stage2_tables(S, bs)
% Stage-2 bar list (22 bars of the two placement-frame sets, frames of the 88-frame bar pattern) and a seeded sequence.
    go = S.g_orient; rows = struct([]);
    labs = {'perpendicular (PD-axis profile)', 'parallel (orthogonal-axis profile)'}; names = {'perp', 'par'};
    for j = 1:2
        g = bs.(['g_' names{j}]); fr = bs.([names{j} '_frames_1based']);
        for rank = 1:numel(fr)
            k = fr(rank);
            row = struct('set', labs{j}, 'bar_orientation_deg', round(go(g + 1), 1), 'bar_pattern_frame_1based', k, 'rank_along_axis', rank, ...
                         'col_center', round(S.cx(k) + 1, 2), 'row_center', round(S.cy(k) + 1, 2));
            if isempty(rows), rows = row; else, rows(end+1) = row; end %#ok<AGROW>
        end
    end
    frames = [rows.bar_pattern_frame_1based]; seq = struct([]); t = 0;
    dur = struct('flash80', [80 920], 'flash14', [14 486]);
    for rep = 1:3
        for cls = {'flash80', 'flash14'}
            d = dur.(cls{1}); perm = frames(randperm(numel(frames)));
            for k = 1:numel(perm)
                row = struct('repeat', rep, 'cls', cls{1}, 'order', k, 'bar_pattern_frame_1based', perm(k), 'flash_ms', d(1), 'gray_ms', d(2), 't_start_s', round(t, 2));
                if isempty(seq), seq = row; else, seq(end+1) = row; end %#ok<AGROW>
                t = t + sum(d) / 1000;
            end
        end
    end
end


%% ===================================================================== geometry (mirror of rapid_prot2_locations.py Frame)
function xy = to_arena(G, s, v)
    p = G.C + s * G.e_pd + v * G.e_or; xy = [p(1) + 1, -p(2) + 1];
end

function sv = to_sv(G, col, row)
    p = [col - 1, -(row - 1)] - G.C; sv = [p * G.e_pd', p * G.e_or'];
end

function q = square(G, s, v, sz)
% Axis-aligned sz-px square whose center is nearest to the PD-frame position (s, v).
    half = (sz - 1) / 2; xy = to_arena(G, s, v); cs = pyround(xy(1) - half); rs = pyround(xy(2) - half);
    sv = to_sv(G, cs + half, rs + half);
    on = G.screen.col_min <= cs && cs + sz - 1 <= G.screen.col_max && G.screen.row_min <= rs && rs + sz - 1 <= G.screen.row_max;
    q = struct('size_px', sz, 'col_start', cs, 'row_start', rs, 'col_center', cs + half, 'row_center', rs + half, ...
               's_px', round(sv(1), 2), 'v_px', round(sv(2), 2), 's_target', s, 'v_target', v, 'on_screen', on);
end

function on = on_screen(G, cs, rs, sz)
    on = G.screen.col_min <= cs && cs + sz - 1 <= G.screen.col_max && G.screen.row_min <= rs && rs + sz - 1 <= G.screen.row_max;
end

function r = pyround(x)
% Python round(): half to even (MATLAB round is half away from zero).
    r = round(x); h = (x - floor(x)) == 0.5; r(h) = 2 * round(x(h) / 2);
end

function [kept, grid] = lattice(G)
    sz = 4; step = 2; half = 1.5; W = struct('s_lead', -12, 's_trail', 16, 'v_half', 8);
    c = to_arena(G, 0, 0); col0 = pyround(c(1) - half); row0 = pyround(c(2) - half);
    reach = ceil(max([abs(W.s_lead), W.s_trail, W.v_half]) * sqrt(2) / step) + 2;
    kept = struct([]); n_off = 0;
    for i = -reach:reach
        for j = -reach:reach
            cs = col0 + i * step; rs = row0 + j * step; sv = to_sv(G, cs + half, rs + half); s = sv(1); v = sv(2);
            if ~(W.s_lead - step / 2 < s && s < W.s_trail + step / 2 && abs(v) < W.v_half + step / 2), continue; end
            q = square(G, s, v, sz);
            if q.col_start ~= cs || q.row_start ~= rs
                q.col_start = cs; q.row_start = rs; q.col_center = cs + half; q.row_center = rs + half; q.s_px = round(s, 2); q.v_px = round(v, 2);
            end
            q.on_screen = on_screen(G, cs, rs, sz);
            if ~q.on_screen, n_off = n_off + 1; continue; end
            q.gen_frame_off = NaN; q.gen_frame_on = NaN;
            if isempty(kept), kept = q; else, kept(end+1) = q; end %#ok<AGROW>
        end
    end
    if isempty(kept), error('rapid_readout:Lattice', 'no lattice positions inside the window'); end
    cols = unique([kept.col_start]); rows = unique([kept.row_start]);
    fl_rows = (rows(end) - rows(1)) / step + 1; fl_cols = (cols(end) - cols(1)) / step + 1; n_fl = fl_rows * fl_cols;
    for k = 1:numel(kept)
        r_i = (kept(k).row_start - rows(1)) / step + 1; c_i = (kept(k).col_start - cols(1)) / step + 1;
        kept(k).gen_frame_off = 1 + (c_i - 1) * fl_rows + r_i; kept(k).gen_frame_on = kept(k).gen_frame_off + n_fl;
    end
    grid = struct('px_rng', struct('row_start', rows(1), 'row_end', rows(end) + sz - 1, 'col_start', cols(1), 'col_end', cols(end) + sz - 1), ...
                  'fl_rows', fl_rows, 'fl_cols', fl_cols, 'n_grid', n_fl, 'n_kept', numel(kept), 'n_dropped_off_screen', n_off, ...
                  'overlap', 1 - step / sz, 'step_px', step, 'size_px', sz);
end

function [placed, notes] = place_sites(G)
    SITES = {0, 0, 'center'; 0, 2, 'center'; -4, 0, 'leading'; -8, 0, 'leading'; 4, 0, 'trailing'; 8, 0, 'trailing'; 12, 0, 'trailing'; ...
             0, 4, 'orthogonal'; 0, -4, 'orthogonal'; 0, 14, 'outside'; 0, -14, 'outside'};
    MIRROR = {[0 14], [0 -14]; [0 -14], [0 14]; [0 4], [0 -4]; [0 -4], [0 4]};
    placed = struct([]); notes = {};
    for i = 1:size(SITES, 1)
        s = SITES{i, 1}; v = SITES{i, 2}; reg = SITES{i, 3};
        q = square(G, s, v, 4); q.region = reg; q.note = '';
        if ~q.on_screen
            m = []; for k = 1:size(MIRROR, 1), if isequal(MIRROR{k, 1}, [s v]), m = MIRROR{k, 2}; end, end
            if ~isempty(m)
                qm = square(G, m(1), m(2), 4); qm.region = reg; qm.note = '';
            else
                qm = struct('on_screen', false);
            end
            if qm.on_screen
                q = qm; q.s_target = s; q.v_target = v;
                q.note = sprintf('mirrored substitute: (%+d, %+d) is off screen, square placed at (%+d, %+d)', s, v, m(1), m(2));
                notes{end+1} = q.note; %#ok<AGROW>
            else
                notes{end+1} = sprintf('site (%+d, %+d) %s is off screen and has no on-screen mirror: dropped', s, v, reg); %#ok<AGROW>
                continue;
            end
        end
        if isempty(placed), placed = q; else, placed(end+1) = q; end %#ok<AGROW>
    end
end

function out = dedupe(sqs)
    out = struct([]); seen = {};
    for k = 1:numel(sqs)
        q = sqs(k); key = sprintf('%d_%d_%d', q.size_px, q.col_start, q.row_start);
        if any(strcmp(seen, key)) || ~q.on_screen, continue; end
        seen{end+1} = key; %#ok<AGROW>
        if isempty(out), out = q; else, out(end+1) = q; end %#ok<AGROW>
    end
end

function o = p3b_order()
    o = repmat([0.5, -0.5, 0, -0.5, 0.5], 1, 4);
end


%% ===================================================================== protocol tables (mirror of rapid_prot2_locations.py)
function [frames, blocks, t] = prot3b_tables(G, sites)
    fld = {'size_px', 'col_start', 'row_start', 'col_center', 'row_center', 's_px', 'v_px', 's_target', 'v_target', 'region'};
    frames = frame_row(1, 'background', 'gray', [], fld);
    frames.size_px = 0;
    pols = {'OFF', 'ON'};
    for p = 1:2
        for k = 1:numel(sites)
            frames(end+1) = frame_row(numel(frames) + 1, 'site', pols{p}, sites(k), fld); %#ok<AGROW>
        end
    end
    hold = [];
    for s = [0, -4, -8]
        q = square(G, s, 0, 3); q.region = 'mi9_series';
        if ~q.on_screen, continue; end
        hold(end+1) = numel(frames) + 1; %#ok<AGROW>
        frames(end+1) = frame_row(numel(frames) + 1, 'mi9_off', 'OFF', q, fld); %#ok<AGROW>
    end
    site_frames = [frames(strcmp({frames.item}, 'site')).frame];
    blocks = struct([]); t = 0; order = p3b_order();
    for b = 1:20
        seq = site_frames(randperm(numel(site_frames))); dur = 3.0 + numel(seq) * 0.6;
        blk = struct('block', b, 'phase', 'rin', 'current_pA', order(b), 'ao_volts', round(order(b) * 0.25, 4), 'settle_s', 3.0, 'n_flashes', numel(seq), ...
                     'flash_ms', 160, 'gray_ms', 440, 'duration_s', round(dur, 1), 't_start_s', round(t, 1), 'frames', strjoin(arrayfun(@(x) sprintf('%d', x), seq, 'UniformOutput', false), ' '));
        if isempty(blocks), blocks = blk; else, blocks(end+1) = blk; end %#ok<AGROW>
        t = t + dur;
    end
    hold_order = [-1 0 1 1 0 -1]; pool = repmat(hold, 1, 6);
    for b = 21:26
        seq = pool(randperm(numel(pool))); dur = 3.0 + numel(seq) * 1.0;
        blocks(end+1) = struct('block', b, 'phase', 'mi9_holding', 'current_pA', hold_order(b - 20), 'ao_volts', round(hold_order(b - 20) * 0.25, 4), 'settle_s', 3.0, 'n_flashes', numel(seq), ...
                     'flash_ms', 500, 'gray_ms', 500, 'duration_s', round(dur, 1), 't_start_s', round(t, 1), 'frames', strjoin(arrayfun(@(x) sprintf('%d', x), seq, 'UniformOutput', false), ' ')); %#ok<AGROW>
        t = t + dur;
    end
end

function f = frame_row(num, item, pol, q, fld)
    f = struct('frame', num, 'item', item, 'size_px', NaN, 'polarity', pol);
    for k = 1:numel(fld)
        if isempty(q), if strcmp(fld{k}, 'region'), f.(fld{k}) = ''; else, f.(fld{k}) = NaN; end
        else, f.(fld{k}) = q.(fld{k}); end
    end
end

function [frames, seq, t, counts, per_item_s] = prot4_tables(G, kept)
    fld = {'size_px', 'col_start', 'row_start', 'col_center', 'row_center', 's_px', 'v_px'};
    frames = frame_row(1, 'background', 'gray', [], fld); frames.size_px = 0;
    items = struct('lattice', [], 'core', [], 'line', [], 'off3', []);
    pols = {'OFF', 'ON'};
    for p = 1:2
        for k = 1:numel(kept)
            frames(end+1) = frame_row(numel(frames) + 1, 'lattice', pols{p}, kept(k), fld); items.lattice(end+1) = frames(end).frame; %#ok<AGROW>
        end
    end
    core = struct([]);
    for s = -6:6, for v = -4:4, q = square(G, s, v, 3); if isempty(core), core = q; else, core(end+1) = q; end, end, end %#ok<AGROW>
    core = dedupe(core);
    for k = 1:numel(core), frames(end+1) = frame_row(numel(frames) + 1, 'core', 'ON', core(k), fld); items.core(end+1) = frames(end).frame; end %#ok<AGROW>
    line = struct([]);
    for s = -8:8, q = square(G, s, 0, 2); if isempty(line), line = q; else, line(end+1) = q; end, end %#ok<AGROW>
    for v = -6:6, line(end+1) = square(G, 0, v, 2); end %#ok<AGROW>
    line = dedupe(line);
    for k = 1:numel(line), frames(end+1) = frame_row(numel(frames) + 1, 'line', 'ON', line(k), fld); items.line(end+1) = frames(end).frame; end %#ok<AGROW>
    off3 = struct([]);
    for s = [0, -4, -8, 4], q = square(G, s, 0, 3); if isempty(off3), off3 = q; else, off3(end+1) = q; end, end %#ok<AGROW>
    off3 = dedupe(off3);
    for k = 1:numel(off3), frames(end+1) = frame_row(numel(frames) + 1, 'off3', 'OFF', off3(k), fld); items.off3(end+1) = frames(end).frame; end %#ok<AGROW>
    dur = struct('lattice', [160 440], 'core', [160 240], 'line', [160 240], 'off3', [160 240]);
    cols = [frames.col_center]; rows = [frames.row_center];
    seq = struct([]); t = 0;
    for p = 1:6
        pool_fr = [items.core, items.line, items.off3]; pool_it = [repmat({'core'}, 1, numel(items.core)), repmat({'line'}, 1, numel(items.line)), repmat({'off3'}, 1, numel(items.off3))];
        if any(p == [1 3 5]), pool_fr = [pool_fr, items.lattice]; pool_it = [pool_it, repmat({'lattice'}, 1, numel(items.lattice))]; end
        perm = randperm(numel(pool_fr)); pool_fr = pool_fr(perm); pool_it = pool_it(perm);
        last = []; k = 0;
        while ~isempty(pool_fr)
            if isempty(last), i = 1;
            else
                ok = find(hypot(cols(pool_fr) - last(1), rows(pool_fr) - last(2)) >= 10, 1);
                if isempty(ok), i = 1; else, i = ok; end
            end
            fr = pool_fr(i); it = pool_it{i}; pool_fr(i) = []; pool_it(i) = [];
            last = [cols(fr), rows(fr)]; k = k + 1; f = frames(fr); d = dur.(it);
            row = struct('pass', p, 'order', k, 'frame', fr, 'item', it, 'size_px', f.size_px, 'polarity', f.polarity, 'col_start', f.col_start, 'row_start', f.row_start, ...
                         's_px', f.s_px, 'v_px', f.v_px, 'flash_ms', d(1), 'gray_ms', d(2), 't_start_s', round(t, 2));
            if isempty(seq), seq = row; else, seq(end+1) = row; end %#ok<AGROW>
            t = t + sum(d) / 1000;
        end
    end
    counts = struct('lattice', numel(items.lattice), 'core', numel(items.core), 'line', numel(items.line), 'off3', numel(items.off3));
    per_item_s = struct('lattice', numel(items.lattice) * 0.6 * 3, 'core', numel(items.core) * 0.4 * 6, 'line', numel(items.line) * 0.4 * 6, 'off3', numel(items.off3) * 0.4 * 6);
end


%% ===================================================================== validation, output helpers
function val = validate_metrics(A, opts, sid)
    val = [];
    if isempty(opts.Metrics), return; end
    T = readtable(opts.Metrics, 'TextType', 'char');
    rows = find(strcmp(T.cell, opts.Cell) | strcmp(T.session, sid));
    if isempty(rows), val = struct('status', sprintf('no row for %s in %s', opts.Cell, opts.Metrics), 'detail', ''); return; end
    m = T(rows(1), :);
    names = {'pd_deg', 'pd_vector_deg', 'm6_pd_1based', 'm6_or_1based', 'g_pd', 'g_or'};
    ours = [A.pd, A.pd_vector, A.m6_pd, A.m6_or, A.g_pd, A.g_or]; tol = [1e-6, 0.051, 0, 0, 0, 0];
    okk = true(1, 6); detail = '';
    for k = 1:6
        ref = m.(names{k});
        if isnan(ours(k)), val.checks.(names{k}) = struct('ours', NaN, 'metrics', ref, 'ok', true, 'note', 'TUTL set not presented'); continue; end
        okk(k) = abs(ours(k) - ref) <= tol(k);
        val.checks.(names{k}) = struct('ours', ours(k), 'metrics', ref, 'ok', okk(k));
        if ~okk(k), detail = [detail sprintf(' %s ours %g metrics %g;', names{k}, ours(k), ref)]; end %#ok<AGROW>
    end
    if all(okk), val.status = 'PASS'; detail = 'all six fields equal'; else, val.status = 'FAIL'; end
    val.detail = strtrim(detail); val.row_cell = m.cell{1};
end

function write_csv(path, rows, cols, int_cols)
% CSV in the Python csv module's style: given column order, True/False booleans, empty for missing, repr-like floats.
    fid = fopen(path, 'w'); fprintf(fid, '%s\n', strjoin(cols, ','));
    for r = 1:numel(rows)
        parts = cell(1, numel(cols));
        for c = 1:numel(cols)
            if isfield(rows, cols{c}), v = rows(r).(cols{c}); else, v = []; end
            parts{c} = fmt_field(v, any(strcmp(int_cols, cols{c})));
        end
        fprintf(fid, '%s\n', strjoin(parts, ','));
    end
    fclose(fid);
end

function s = fmt_field(v, is_int)
    if ischar(v) || isstring(v)
        s = char(v);
        if any(s == ',') || any(s == '"') || any(s == newline), s = ['"' strrep(s, '"', '""') '"']; end
    elseif islogical(v)
        if v, s = 'True'; else, s = 'False'; end
    elseif isempty(v) || (isnumeric(v) && isscalar(v) && isnan(v))
        s = '';
    elseif is_int
        s = sprintf('%d', v);
    elseif v == round(v) && abs(v) < 1e15
        s = sprintf('%.1f', v);
    else
        s = sprintf('%.15g', v);
        if str2double(s) ~= v, s = sprintf('%.17g', v); end
    end
end

function c = crc32(str)
% CRC-32 (IEEE 802.3, as zlib.crc32) of a char vector, pure MATLAB.
    b = uint8(str); c = uint32(4294967295);
    for i = 1:numel(b)
        c = bitxor(c, uint32(b(i)));
        for j = 1:8
            if bitand(c, uint32(1)), c = bitxor(bitshift(c, -1), uint32(3988292384)); else, c = bitshift(c, -1); end
        end
    end
    c = double(bitxor(c, uint32(4294967295)));
end
