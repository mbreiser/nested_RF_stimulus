classdef SessionProgress < handle
    %SESSIONPROGRESS Where run_session is: the step, the row, the time left, in the command window and a small window.
    %   P = SessionProgress(S, show_window) opens the status window when show_window is true and MATLAB has a
    %   display; the text lines are printed either way. run_session calls
    %     P.step_start(S, k)                         a step begins (prints the step and session time left)
    %     P.row(k, n, label, done_s, total_s)        before row k of n of the current step (returns a text suffix)
    %     P.countdown(left_s, total_s, label)        a timed wait without rows (the dark period)
    %     P.step_end(S, k)                           a step ended (any status)
    %     P.refresh(S)                               the step table changed (estimates, status)
    %   Time left in a step: the programmed time of the rows still to run, scaled by the measured
    %   elapsed / programmed ratio of the rows already run (the controller overhead per row), clamped to
    %   1..1.5. Time left in the session: that plus the estimates of the selected steps still pending.
    %   The window is updated between rows only: during a row MATLAB waits on the G4 Host.
    %   27 October 2026.

    properties
        fig = [];                 % uifigure, empty when there is no window
        h = struct();             % its components
        S = [];                   % last session state seen
        k = 0;                    % current step index
        t_step = [];              % tic of the current step
        step_left_s = 0;          % estimated time left in the current step
        selected = {};            % step ids selected for this run
    end

    methods
        function self = SessionProgress(S, show_window, selected)
            self.S = S;
            if nargin < 3, selected = {S.steps.id}; end
            self.selected = cellstr(selected);
            if show_window && usejava('jvm') && feature('ShowFigureWindows')   % a display (desktop or -nodesktop), not -batch
                try
                    self.build();
                    self.refresh(S);
                catch err
                    fprintf('    (status window not available: %s)\n', err.message);
                    self.fig = [];
                end
            end
        end

        function step_start(self, S, k)
            self.S = S; self.k = k; self.t_step = tic;
            st = S.steps(k);
            self.step_left_s = st.est_min * 60;
            [left_s, t_end] = self.session_left(S, k, self.step_left_s);
            fprintf('    expected %.1f min, ends ~%s | session: %.1f min elapsed, ~%.0f min left, finish ~%s\n', ...
                st.est_min, clock_at(self.step_left_s), self.session_elapsed(S, k) / 60, left_s / 60, t_end);
            if self.has_window()
                self.h.step.Text = sprintf('%d/%d  %s  -  %s', k, numel(S.steps), st.id, st.label);
                self.h.row.Text = sprintf('expected %.1f min, ends ~%s', st.est_min, clock_at(self.step_left_s));
                self.set_bar(0);
                self.refresh(S);
            end
        end

        function txt = row(self, r, n, label, done_s, total_s)
            % before row r: done_s programmed seconds already run in this step, total_s programmed in all
            el = toc(self.t_step);
            if done_s > 0, ratio = min(max(el / done_s, 1), 1.5); else, ratio = 1; end
            self.step_left_s = (total_s - done_s) * ratio;
            [left_s, t_end] = self.session_left(self.S, self.k, self.step_left_s);
            txt = sprintf('| %4.1f/%.1f min, ~%.1f left', el / 60, (el + self.step_left_s) / 60, self.step_left_s / 60);
            if self.has_window()
                self.h.row.Text = sprintf('row %d/%d  %s   |   %.1f of ~%.1f min, ~%.1f min left, step ends ~%s', ...
                    r, n, label, el / 60, (el + self.step_left_s) / 60, self.step_left_s / 60, clock_at(self.step_left_s));
                self.set_bar(done_s / max(total_s, eps));       % programmed time done: right also when the clock runs faster (dry run)
                self.set_session_line(left_s, t_end);
                drawnow limitrate;
            end
        end

        function countdown(self, left_s, total_s, label)
            self.step_left_s = left_s;
            if self.has_window()
                self.h.row.Text = sprintf('%s: %.0f s left', label, left_s);
                self.set_bar(1 - left_s / max(total_s, eps));
                [l, t_end] = self.session_left(self.S, self.k, left_s); self.set_session_line(l, t_end);
                drawnow limitrate;
            end
        end

        function step_end(self, S, k)
            self.S = S; self.step_left_s = 0;
            if self.has_window()
                self.set_bar(1);
                self.h.row.Text = sprintf('%s %s in %.1f min', S.steps(k).id, S.steps(k).status, S.steps(k).elapsed_s / 60);
                self.refresh(S);
            end
        end

        function note(self, msg)
            % the last flag or warning, shown until the next one
            if self.has_window(), self.h.flag.Text = msg; drawnow limitrate; end
        end

        function refresh(self, S)
            self.S = S;
            if ~self.has_window(), return; end
            n = numel(S.steps); D = cell(n, 8);
            ahead = 0;
            for j = 1:n
                st = S.steps(j);
                D{j, 1} = j; D{j, 2} = st.id; D{j, 3} = st.label; D{j, 4} = st.status; D{j, 5} = sprintf('%.1f', st.est_min);
                D{j, 6} = ''; if st.elapsed_s > 0, D{j, 6} = sprintf('%.1f', st.elapsed_s / 60); end
                D{j, 7} = ''; if ~isempty(st.t_start), D{j, 7} = st.t_start(12:min(16, end)); end
                D{j, 8} = '';
                if j == self.k && strcmp(st.status, 'running')
                    ahead = self.step_left_s; D{j, 8} = ['~' clock_at(ahead)];
                elseif strcmp(st.status, 'pending') && ismember(st.id, self.selected)
                    ahead = ahead + st.est_min * 60; D{j, 8} = ['~' clock_at(ahead)];
                elseif ~isempty(st.t_end)
                    D{j, 8} = st.t_end(12:min(16, end));
                end
            end
            self.h.table.Data = D;
            rs = find(strcmp({S.steps.status}, 'running'));
            removeStyle(self.h.table);
            if ~isempty(rs), addStyle(self.h.table, uistyle('BackgroundColor', [1 0.95 0.75], 'FontWeight', 'bold'), 'row', rs); end
            bad = find(ismember({S.steps.status}, {'failed', 'interrupted'}));
            if ~isempty(bad), addStyle(self.h.table, uistyle('FontColor', [0.75 0 0]), 'row', bad); end
            dn = find(strcmp({S.steps.status}, 'done'));
            if ~isempty(dn), addStyle(self.h.table, uistyle('FontColor', [0.35 0.35 0.35]), 'row', dn); end
            [l, t_end] = self.session_left(S, self.k, self.step_left_s); self.set_session_line(l, t_end);
            drawnow limitrate;
        end

        function close(self)
            if self.has_window(), delete(self.fig); end
            self.fig = [];
        end

        function [left_s, t_end] = session_left(self, S, k, step_left_s)
            left_s = step_left_s;
            for j = 1:numel(S.steps)
                if j ~= k && strcmp(S.steps(j).status, 'pending') && ismember(S.steps(j).id, self.selected)
                    left_s = left_s + S.steps(j).est_min * 60;
                end
            end
            t_end = clock_at(left_s);
        end

        function s = session_elapsed(self, S, k)
            % steps already finished plus the current one so far
            s = sum([S.steps([1:k-1, k+1:end]).elapsed_s]);
            if ~isempty(self.t_step), s = s + toc(self.t_step); end
        end
    end

    methods (Access = private)
        function tf = has_window(self)
            tf = ~isempty(self.fig) && isvalid(self.fig);
        end

        function build(self)
            name = sprintf('Session %s', self.S.cell);
            delete(findall(groot, 'Type', 'figure', 'Name', name));   % a resumed or repeated session: one window
            f = uifigure('Name', name, 'Position', [40 80 720 440], 'Color', 'w');
            g = uigridlayout(f, [6 1], 'RowHeight', {24, 22, 18, 22, '1x', 22}, 'Padding', [10 8 10 8], 'RowSpacing', 4);
            self.h.session = uilabel(g, 'Text', sprintf('Session %s%s', self.S.cell, tern(self.S.dry_run, '  [DRY RUN]', '')), ...
                'FontSize', 14, 'FontWeight', 'bold');
            self.h.step = uilabel(g, 'Text', '', 'FontSize', 13);
            bar = uipanel(g, 'BackgroundColor', [0.9 0.9 0.9], 'BorderType', 'none');
            self.h.fill = uipanel(bar, 'BackgroundColor', [0.25 0.55 0.85], 'BorderType', 'none', 'Position', [0 0 1 18]);
            self.h.bar = bar;
            self.h.row = uilabel(g, 'Text', '', 'FontSize', 12);
            self.h.table = uitable(g, 'ColumnName', {'#', 'step', 'what', 'status', 'est min', 'took', 'start', 'end'}, ...
                'ColumnWidth', {26, 62, '1x', 76, 56, 60, 46, 54}, 'RowName', {}, 'FontSize', 12);
            self.h.flag = uilabel(g, 'Text', '', 'FontColor', [0.75 0 0], 'FontSize', 12);
            self.fig = f;
        end

        function set_bar(self, frac)
            frac = min(max(frac, 0), 1);
            w = self.h.bar.InnerPosition(3);
            self.h.fill.Position = [0 0 max(1, round(frac * w)) self.h.bar.InnerPosition(4)];
        end

        function set_session_line(self, left_s, t_end)
            if ~self.has_window(), return; end
            el = self.session_elapsed(self.S, self.k);
            self.h.session.Text = sprintf('Session %s%s  |  %.1f min elapsed, ~%.0f min left, finish ~%s', self.S.cell, ...
                tern(self.S.dry_run, ' [DRY RUN]', ''), el / 60, left_s / 60, t_end);
        end
    end
end

function s = clock_at(seconds_from_now)
s = char(datetime('now') + seconds(seconds_from_now), 'HH:mm');
end

function v = tern(c, a, b)
if c, v = a; else, v = b; end
end
