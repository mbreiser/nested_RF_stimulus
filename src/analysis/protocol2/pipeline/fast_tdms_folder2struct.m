function [Log, info] = fast_tdms_folder2struct(exp_folder, varargin)
% FAST_TDMS_FOLDER2STRUCT  Read a G4 TDMS log folder into the G4_TDMS_folder2struct Log struct, fast.
%
%   LOG = FAST_TDMS_FOLDER2STRUCT(FOLDER) returns the same struct that
%   G4_TDMS_folder2struct(FOLDER) saves as G4_TDMS_Logs_<name>.mat
%   (Log.ADC, Log.AO, Log.Frames, Log.Commands with the same fields,
%   classes, row order and shapes), but in well under a second instead of
%   about a minute, and without the parallel computing toolbox, the NI DLL
%   or .NET.
%
%   [LOG, INFO] = FAST_TDMS_FOLDER2STRUCT(FOLDER, 'Save', true) also writes
%   the .mat where the converter writes it (parent of the TDMS folder,
%   'G4_TDMS_Logs_<folder name>.mat'). INFO holds the resolved TDMS folder,
%   the list of files, the read time per file and the detected channel
%   layout.
%
%   FOLDER may be the session folder (…/<sid>, holding 'Log Files'), the
%   'Log Files' folder (holding one time-stamped subfolder) or the
%   time-stamped TDMS folder itself, as for the converter.
%
%   Why it is fast. The G4 host writes every ADC channel as ONE TDMS segment
%   whose raw-data index says "one value per chunk", so a 15-min session is
%   about 9 million chunks per file. The generic reader
%   (TDMS_readChannelOrGroup) walks the chunks; this function parses the
%   segment lead-in and metadata once and reads the segment's raw data with a
%   single fread, which is disk-bound (about 0.1 s for a 107 MB file from an
%   SSD). The same approach, in Python, is tools/rapid_prot2_fast.py of the
%   ephys-standard-analysis repository (26 Oct 2026), which this port
%   follows; see docs/notes/MATLAB_FAST_READOUT_20261026.md there for the
%   validation record (field-by-field equality with the converter's .mat on
%   25 TUTL sessions of the historic layout and on five October 2026
%   sessions of Jin Yong's rig, both layouts).
%
%   Supported TDMS features: little-endian files; segments with and without
%   metadata, carried-over object lists, "same as previous" raw-data indices;
%   non-interleaved and interleaved numeric raw data; string channels (the
%   command log); several segments per file. Not supported (error): big-endian
%   files, DAQmx raw data.
%
%   Channel layouts. The function is layout-agnostic like the converter:
%   Log.ADC rows follow the *_ADC<n>_* files in name order (historic rig:
%   ADC0 = frame index, ADC1 = Vm; Jin Yong's rig from October 2026: ADC1 =
%   Vm, ADC2 = current monitor, with the frame index in Log.Frames).
%   INFO.layout names the layout from the data (first ADC channel integer-
%   valued => historic).
%
%   Options (name/value): 'Save' (false), 'Verbose' (false).
%
%   See also G4_TDMS_FOLDER2STRUCT, RAPID_READOUT.

    opts = struct('Save', false, 'Verbose', false);
    for k = 1:2:numel(varargin)
        opts.(varargin{k}) = varargin{k+1};
    end
    t_all = tic;

    %% Resolve the TDMS folder (session folder / 'Log Files' / the folder itself)
    tdms_folder = resolve_tdms_folder(exp_folder);
    parent = fileparts(tdms_folder);
    [~, folder_name] = fileparts(tdms_folder);

    files = dir(fullfile(tdms_folder, '*.tdms'));
    files = files(~[files.isdir]);
    names = sort({files.name});             % the converter's dir() order
    assert(~isempty(names), 'fast_tdms_folder2struct:NoFiles', 'cannot find any .tdms files in %s', tdms_folder);

    %% Read every file (same classification as G4_TDMS_folder2struct)
    n = numel(names);
    ADCTime = {}; ADCVolts = {}; ADCChannels = {};
    AOTime = {}; AOVolts = {}; AOChannels = {};
    FrameTime = {}; FramePosition = {};
    CommandTime = []; CommandName = {}; CommandData = {};
    info = struct('tdms_folder', tdms_folder, 'files', {names}, 'read_seconds', zeros(1, n), 'layout', '', 'total_seconds', NaN, 'empty_files', {{}});

    for i = 1:n
        fn = names{i};
        path = fullfile(tdms_folder, fn);
        t0 = tic;
        if contains(fn, 'ADC')
            chan = str2double(fn(strfind(fn, 'ADC') + 3));
            group = ['ADC' num2str(chan)];
            [chans, data] = read_tdms_group(path, group);
            x = first_data(chans, data, path, group);
            if contains(fn, 'Time')
                ADCTime{end+1} = x(:)';  %#ok<AGROW>
                ADCChannels{end+1} = group; %#ok<AGROW>
            elseif contains(fn, 'Volts')
                ADCVolts{end+1} = x(:)'; %#ok<AGROW>
            else
                error('fast_tdms_folder2struct:Unexpected', 'unexpected TDMS file: neither volts nor time: %s', fn);
            end
        elseif contains(fn, 'AO')
            chan = str2double(fn(strfind(fn, 'AO') + 2));
            group = ['AO' num2str(chan - 2)];
            [chans, data] = read_tdms_group(path, group);
            if isempty(chans)
                % an AO channel set active but never driven leaves a header-only file (Prot3, 2 Oct 2026): no AO data
                info.empty_files{end+1} = fn; info.read_seconds(i) = toc(t0);
                continue
            end
            x = first_data(chans, data, path, group);
            if contains(fn, 'Time')
                AOTime{end+1} = x(:)'; %#ok<AGROW>
                AOChannels{end+1} = group; %#ok<AGROW>
            elseif contains(fn, 'Volts')
                AOVolts{end+1} = x(:)'; %#ok<AGROW>
            else
                error('fast_tdms_folder2struct:Unexpected', 'unexpected TDMS file: neither volts nor time: %s', fn);
            end
        elseif contains(fn, 'Frame')
            group = 'Pattern Position';
            [chans, data] = read_tdms_group(path, group);
            if isempty(chans)
                % no pattern displayed while logging (the dark period, allOff) leaves header-only frame files
                info.empty_files{end+1} = fn; info.read_seconds(i) = toc(t0);
                continue
            end
            x = first_data(chans, data, path, group);   % 'X Index' (position) or 'Time'
            if contains(fn, 'Position')
                FramePosition{end+1} = x(:)'; %#ok<AGROW>
            elseif contains(fn, 'Time')
                FrameTime{end+1} = x(:)'; %#ok<AGROW>
            else
                error('fast_tdms_folder2struct:Unexpected', 'unexpected TDMS file: neither pattern position nor time: %s', fn);
            end
        else
            group = 'Commands Received';
            [chans, data] = read_tdms_group(path, group);
            if isempty(chans)
                error('fast_tdms_folder2struct:Unexpected', ['Unexpected file in TDMS folder: ' fn '. Could not find the "Commands Received" group in this file']);
            end
            CommandTime = pick(chans, data, 'Time'); CommandTime = CommandTime(:)';
            CommandName = pick(chans, data, 'Command'); CommandName = CommandName(:)';
            CommandData = pick(chans, data, 'Data'); CommandData = CommandData(:)';
        end
        info.read_seconds(i) = toc(t0);
        if opts.Verbose, fprintf('  %-45s %6.3f s\n', fn, info.read_seconds(i)); end
    end

    %% Assemble the struct exactly as the converter does
    missing = '';
    if isempty(ADCTime), missing = [missing ' ADCTime']; end
    if isempty(ADCVolts), missing = [missing ' ADCVolts']; end
    if isempty(ADCChannels), missing = [missing ' ADCChannels']; end
    if isempty(AOTime), missing = [missing ' AOTime']; end
    if isempty(AOVolts), missing = [missing ' AOVolts']; end
    if isempty(AOChannels), missing = [missing ' AOChannels']; end
    if isempty(FrameTime), missing = [missing ' FrameTime']; end
    if isempty(FramePosition), missing = [missing ' FramePosition']; end
    if isempty(CommandTime), missing = [missing ' CommandTime']; end
    if isempty(CommandName), missing = [missing ' CommandName']; end
    if isempty(CommandData), missing = [missing ' CommandData']; end
    if ~isempty(missing) && opts.Verbose
        fprintf(['variables not found/imported: ' missing '\n']);
    end

    Log.ADC.Time = stack_rows(ADCTime);
    Log.ADC.Volts = stack_rows(ADCVolts);
    Log.ADC.Channels = cell_or_empty(ADCChannels);
    Log.AO.Time = stack_rows(AOTime);
    Log.AO.Volts = stack_rows(AOVolts);
    Log.AO.Channels = cell_or_empty(AOChannels);
    Log.Frames.Time = stack_rows(FrameTime);
    Log.Frames.Position = stack_rows(FramePosition);
    Log.Commands.Time = CommandTime;
    Log.Commands.Name = CommandName;
    Log.Commands.Data = CommandData;

    % Layout from the data, as pharm_protocol.load_session decides it
    if ~isempty(Log.ADC.Volts)
        v0 = Log.ADC.Volts(1, :);
        if mean(abs(v0 - round(v0)) < 1e-3) > 0.9
            info.layout = 'historic: ADC0 frame index / ADC1 Vm';
        else
            info.layout = 'October 2026: ADC1 Vm / ADC2 current + frame log';
        end
    end
    info.total_seconds = toc(t_all);

    if opts.Save
        out_path = fullfile(parent, ['G4_TDMS_Logs_' folder_name '.mat']);
        save(out_path, 'Log');
        info.saved_to = out_path;
        if opts.Verbose, fprintf('saved %s\n', out_path); end
    end
end


%% ============================ folder and struct helpers ============================

function tdms_folder = resolve_tdms_folder(exp_folder)
    exp_folder = char(exp_folder);
    if ~isempty(exp_folder) && any(exp_folder(end) == '/\')
        exp_folder = exp_folder(1:end-1);
    end
    assert(isfolder(exp_folder), 'fast_tdms_folder2struct:NoFolder', 'not a folder: %s', exp_folder);
    if ~isempty(dir(fullfile(exp_folder, '*.tdms')))
        tdms_folder = exp_folder; return;
    end
    if isfolder(fullfile(exp_folder, 'Log Files'))
        tdms_folder = resolve_tdms_folder(fullfile(exp_folder, 'Log Files')); return;
    end
    d = dir(exp_folder);
    d = d([d.isdir] & ~ismember({d.name}, {'.', '..'}));
    for k = 1:numel(d)
        if ~isempty(dir(fullfile(exp_folder, d(k).name, '*.tdms')))
            tdms_folder = fullfile(exp_folder, d(k).name); return;
        end
    end
    error('fast_tdms_folder2struct:NoFiles', 'no TDMS log folder under %s', exp_folder);
end

function x = first_data(chans, data, path, group)
    assert(~isempty(chans), 'fast_tdms_folder2struct:NoGroup', 'group "%s" has no data in %s', group, path);
    x = data{1};
end

function x = pick(chans, data, name)
    k = find(strcmp(chans, name), 1);
    assert(~isempty(k), 'fast_tdms_folder2struct:NoChannel', 'channel "%s" not found', name);
    x = data{k};
end

function M = stack_rows(rows)
    if isempty(rows)
        M = [];
    else
        M = vertcat(rows{:});
    end
end

function c = cell_or_empty(c)
    if isempty(c), c = zeros(1, 0); end   % matches the converter: [](logical(zeros(1,n))) -> 1x0 double
end


%% ================================ TDMS reader ================================

function [chan_names, chan_data] = read_tdms_group(path, group)
% Read every channel of GROUP that carries raw data. CHAN_DATA{k} is a column
% vector (numeric class as stored, or a cell array of char for string channels),
% with the segments concatenated.
    TOC_META = 2; TOC_NEWOBJ = 4; TOC_RAW = 8; TOC_INTERLEAVED = 32; TOC_BIGENDIAN = 64; TOC_DAQMX = 128;
    STRING = 32;

    fid = fopen(path, 'r', 'ieee-le');
    assert(fid > 0, 'fast_tdms_folder2struct:Open', 'cannot open %s', path);
    cleaner = onCleanup(@() fclose(fid));
    fseek(fid, 0, 'eof'); fsize = ftell(fid); fseek(fid, 0, 'bof');

    cur = cell(0, 2);         % current object list: {path, raw} per row ({'/group/chan', [dtype nvals totalbytes] or []})
    raw_prev = containers.Map('KeyType', 'char', 'ValueType', 'any');
    out = containers.Map('KeyType', 'char', 'ValueType', 'any');   % key -> cell of parts
    order = {};
    off = 0;
    while off + 28 <= fsize
        fseek(fid, off, 'bof');
        tag = fread(fid, 4, 'uint8=>char')';
        toc = fread(fid, 1, 'uint32=>double');
        fread(fid, 1, 'uint32=>double');                       % version
        nextseg = fread(fid, 1, 'uint64=>uint64');
        rawoff = fread(fid, 1, 'uint64=>double');
        if ~strcmp(tag, 'TDSm')
            error('fast_tdms_folder2struct:Format', 'not a TDMS segment at byte %d of %s', off, path);
        end
        if bitand(toc, TOC_BIGENDIAN), error('fast_tdms_folder2struct:Format', 'big-endian TDMS is not supported (%s)', path); end
        if bitand(toc, TOC_DAQMX), error('fast_tdms_folder2struct:Format', 'DAQmx TDMS is not supported (%s)', path); end
        if nextseg == intmax('uint64')
            seg_end = fsize;                                   % incomplete last segment: to end of file
        else
            seg_end = min(off + 28 + double(nextseg), fsize);
        end
        data_start = off + 28 + rawoff;

        if bitand(toc, TOC_META)
            meta = fread(fid, rawoff, 'uint8=>uint8');
            objs = parse_meta(meta);                            % n x 3 cell: path, raw ([] / 'prev' / [dtype nvals tb]), props
            if bitand(toc, TOC_NEWOBJ), cur = cell(0, 2); end
            for k = 1:size(objs, 1)
                pth = objs{k, 1}; raw = objs{k, 2};
                if ischar(raw) && strcmp(raw, 'prev')
                    if isKey(raw_prev, pth), raw = raw_prev(pth); else, raw = []; end
                elseif ~isempty(raw)
                    raw_prev(pth) = raw;
                end
                hit = find(strcmp(cur(:, 1), pth), 1);
                if isempty(hit)
                    cur(end+1, :) = {pth, raw}; %#ok<AGROW>
                else
                    cur{hit, 2} = raw;
                end
            end
        end
        if isempty(cur)
            chans = cell(0, 2);
        else
            has = ~cellfun(@isempty, cur(:, 2));
            chans = cur(has, :);
        end
        if bitand(toc, TOC_RAW) && ~isempty(chans) && seg_end > data_start
            nbytes = seg_end - data_start;
            fseek(fid, data_start, 'bof');
            nch = size(chans, 1);
            dtypes = cellfun(@(r) r(1), chans(:, 2));
            nvals = cellfun(@(r) r(2), chans(:, 2));
            if any(dtypes == STRING)
                % one chunk = per channel [uint32 end offsets x nvals][string bytes] (total bytes given) or fixed numeric records
                chunk_bytes = 0;
                for k = 1:nch
                    r = chans{k, 2};
                    if r(1) == STRING, chunk_bytes = chunk_bytes + r(3); else, chunk_bytes = chunk_bytes + item_size(r(1)) * r(2); end
                end
                nchunks = floor(nbytes / chunk_bytes);
                buf = fread(fid, nchunks * chunk_bytes, 'uint8=>uint8');
                pos = 0;
                for c = 1:nchunks
                    for k = 1:nch
                        r = chans{k, 2}; key = chans{k, 1};
                        if r(1) == STRING
                            ends = double(typecast(buf(pos + 1:pos + 4 * r(2)), 'uint32')); pos = pos + 4 * r(2);
                            nb = r(3) - 4 * r(2);
                            s = buf(pos + 1:pos + nb); pos = pos + nb;
                            starts = [0; ends(1:end-1)];
                            vals = cell(r(2), 1);
                            for j = 1:r(2)
                                vals{j} = bytes_to_char(s(starts(j) + 1:ends(j)));
                            end
                        else
                            [mtype, isz] = mat_type(r(1));
                            vals = typecast(buf(pos + 1:pos + isz * r(2)), mtype); vals = vals(:); pos = pos + isz * r(2);
                        end
                        [out, order] = append_part(out, order, key, vals);
                    end
                end
            elseif bitand(toc, TOC_INTERLEAVED)
                if numel(unique(dtypes)) == 1 && all(nvals == 1)
                    [mtype, isz] = mat_type(dtypes(1));
                    nrec = floor(nbytes / (isz * nch));
                    arr = fread(fid, [nch, nrec], [mtype '=>' mtype]);
                    for k = 1:nch
                        [out, order] = append_part(out, order, chans{k, 1}, arr(k, :)');
                    end
                else
                    [out, order] = read_records(fid, nbytes, chans, out, order);
                end
            else
                if nch == 1
                    [mtype, isz] = mat_type(dtypes(1));
                    vals = fread(fid, floor(nbytes / isz), [mtype '=>' mtype]);
                    [out, order] = append_part(out, order, chans{1, 1}, vals(:));
                else
                    [out, order] = read_records(fid, nbytes, chans, out, order);
                end
            end
        end
        off = seg_end;
    end

    % Assemble the requested group, channels in file order
    chan_names = {}; chan_data = {};
    prefix = ['/''' strrep(group, '''', '''''') '''/'''];
    for k = 1:numel(order)
        key = order{k};
        if strncmp(key, prefix, numel(prefix))
            cname = strrep(key(numel(prefix) + 1:end-1), '''''', '''');
            parts = out(key);
            v = vertcat(parts{:});                 % numeric parts, or cells of char for string channels
            chan_names{end+1} = cname; %#ok<AGROW>
            chan_data{end+1} = v; %#ok<AGROW>
        end
    end
end

function [out, order] = read_records(fid, nbytes, chans, out, order)
% Non-interleaved chunk with several channels (or interleaved-by-value with unequal counts): fixed records of
% [ch1 x nvals1][ch2 x nvals2]... read as bytes and sliced.
    nch = size(chans, 1);
    sizes = zeros(nch, 1); types = cell(nch, 1);
    for k = 1:nch
        r = chans{k, 2}; [types{k}, isz] = mat_type(r(1)); sizes(k) = isz * r(2);
    end
    rec = sum(sizes); nrec = floor(nbytes / rec);
    buf = fread(fid, [rec, nrec], 'uint8=>uint8');
    p = 0;
    for k = 1:nch
        b = buf(p + 1:p + sizes(k), :); p = p + sizes(k);
        vals = typecast(b(:), types{k});
        [out, order] = append_part(out, order, chans{k, 1}, vals(:));
    end
end

function [out, order] = append_part(out, order, key, vals)
    if isKey(out, key)
        parts = out(key); parts{end+1} = vals; out(key) = parts;
    else
        out(key) = {vals}; order{end+1} = key;
    end
end

function objs = parse_meta(buf)
% Parse one segment's metadata: {path, raw, props} per object. raw = [] (no raw data), 'prev' (same as the
% previous segment) or [dtype nvals total_bytes] (total_bytes NaN unless string).
    STRING = 32; TIMESTAMP = 68;
    p = 0;
    nobj = double(typecast(buf(p + 1:p + 4), 'uint32')); p = p + 4;
    objs = cell(nobj, 3);
    for k = 1:nobj
        L = double(typecast(buf(p + 1:p + 4), 'uint32')); p = p + 4;
        pth = bytes_to_char(buf(p + 1:p + L)); p = p + L;
        idx = double(typecast(buf(p + 1:p + 4), 'uint32')); p = p + 4;
        if idx == 0
            raw = 'prev';
        elseif idx == 4294967295
            raw = [];
        else
            if any(idx == [hex2dec('69120000') hex2dec('69130000') hex2dec('1269') hex2dec('1369')])
                error('fast_tdms_folder2struct:Format', 'DAQmx raw data is not supported');
            end
            dt = double(typecast(buf(p + 1:p + 4), 'uint32')); p = p + 4;
            p = p + 4;                                                    % array dimension (always 1)
            nv = double(typecast(buf(p + 1:p + 8), 'uint64')); p = p + 8;
            tb = NaN;
            if dt == STRING, tb = double(typecast(buf(p + 1:p + 8), 'uint64')); p = p + 8; end
            raw = [dt nv tb];
        end
        nprop = double(typecast(buf(p + 1:p + 4), 'uint32')); p = p + 4;
        props = struct();
        for j = 1:nprop
            L = double(typecast(buf(p + 1:p + 4), 'uint32')); p = p + 4;
            name = bytes_to_char(buf(p + 1:p + L)); p = p + L;
            t = double(typecast(buf(p + 1:p + 4), 'uint32')); p = p + 4;
            if t == STRING
                L = double(typecast(buf(p + 1:p + 4), 'uint32')); p = p + 4;
                val = bytes_to_char(buf(p + 1:p + L)); p = p + L;
            elseif t == TIMESTAMP
                frac = typecast(buf(p + 1:p + 8), 'uint64'); sec = typecast(buf(p + 9:p + 16), 'int64'); p = p + 16;
                val = [double(sec) double(frac)];
            else
                [mtype, isz] = mat_type(t);
                val = typecast(buf(p + 1:p + isz), mtype); p = p + isz;
            end
            props.(matlab.lang.makeValidName(name)) = val;
        end
        objs(k, :) = {pth, raw, props};
    end
end

function [mtype, isz] = mat_type(code)
    switch code
        case 1,  mtype = 'int8';   isz = 1;
        case 2,  mtype = 'int16';  isz = 2;
        case 3,  mtype = 'int32';  isz = 4;
        case 4,  mtype = 'int64';  isz = 8;
        case 5,  mtype = 'uint8';  isz = 1;
        case 6,  mtype = 'uint16'; isz = 2;
        case 7,  mtype = 'uint32'; isz = 4;
        case 8,  mtype = 'uint64'; isz = 8;
        case 9,  mtype = 'single'; isz = 4;
        case 10, mtype = 'double'; isz = 8;
        case 33, mtype = 'uint8';  isz = 1;      % boolean
        otherwise, error('fast_tdms_folder2struct:Format', 'unsupported TDMS data type code %d', code);
    end
end

function isz = item_size(code)
    [~, isz] = mat_type(code);
end

function s = bytes_to_char(b)
    if isempty(b)
        s = '';
    else
        s = native2unicode(reshape(uint8(b), 1, []), 'UTF-8');
    end
end
