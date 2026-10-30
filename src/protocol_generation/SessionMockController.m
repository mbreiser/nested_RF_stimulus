classdef SessionMockController < handle
    %SESSIONMOCKCONTROLLER Offline stand-in for PanelsController used by run_session's dry run.
    %   No TCP, no G4 Host, no arena. Method names and argument validators follow
    %   PanelsController (G4_Display_Tools PControl_Matlab/controller) so that a
    %   sequence that passes here is at least well-formed for the real controller.
    %   Every call is recorded in `calls` ({method, args}); startDisplay waits
    %   `sim_speed` x the programmed duration (default 0: no waiting).
    %
    %   Failure injection for tests: set fail_method / fail_on_call (nth call of
    %   that method) / fail_mode 'error' | 'false'.
    %
    %   Derived from Prot3MockController (TUTL p4_experiment_design, 2 Oct 2026).
    %   26 October 2026.

    properties
        calls = {};
        isOpen = false;
        logRunning = false;
        root = '';
        sim_speed = 0;          % fraction of real time to wait in startDisplay
        verbose = false;
        fail_method = '';
        fail_on_call = 1;
        fail_mode = 'error';
        counts = struct();
    end

    methods
        function self = SessionMockController(opts)
            if nargin > 0 && isstruct(opts)
                for f = fieldnames(opts).', self.(f{1}) = opts.(f{1}); end
            end
        end

        function open(self, startHost)
            arguments
                self (1,1) SessionMockController
                startHost (1,1) logical = false
            end
            self.rec('open', {startHost});
            self.isOpen = true;
        end
        function close(self, stopHost)
            arguments
                self (1,1) SessionMockController
                stopHost (1,1) logical = false
            end
            self.rec('close', {stopHost});
            self.isOpen = false;
        end
        function rtn = stopDisplay(self)
            rtn = self.rec('stopDisplay', {});
        end
        function rtn = allOff(self)
            rtn = self.rec('allOff', {});
        end
        function rtn = allOn(self)
            rtn = self.rec('allOn', {});
        end
        function rtn = setRootDirectory(self, dirName, createDir)
            arguments
                self (1,1) SessionMockController
                dirName (1,1) string
                createDir (1,1) logical = true
            end
            self.root = char(dirName);
            rtn = self.rec('setRootDirectory', {char(dirName), createDir});
        end
        function rtn = setActiveAOChannels(self, ch)
            arguments
                self (1,1) SessionMockController
                ch (1,:) {mustBeInteger, mustBeGreaterThanOrEqual(ch, 2), mustBeLessThanOrEqual(ch, 5)} = [2 3 4 5]
            end
            rtn = self.rec('setActiveAOChannels', {ch});
        end
        function rtn = setActiveAIChannels(self, ch)
            arguments
                self (1,1) SessionMockController
                ch (1,:) {mustBeInteger, mustBeGreaterThanOrEqual(ch, 0), mustBeLessThanOrEqual(ch, 3)} = [0 1 2 3]
            end
            rtn = self.rec('setActiveAIChannels', {ch});
        end
        function rtn = startLog(self)
            rtn = self.rec('startLog', {});
            if rtn, self.logRunning = true; end
        end
        function rtn = stopLog(self, options)
            arguments
                self (1,1) SessionMockController
                options.timeout (1,1) double {mustBeGreaterThan(options.timeout, 0)} = 30.0
                options.showTimeoutMessage (1,1) logical = false
                options.showTimeoutDialog (1,1) logical = false
                options.timeoutMessage (1,1) string = "Log failed to stop. Please stop manually."
            end
            rtn = self.rec('stopLog', {options.showTimeoutMessage});
            if rtn, self.logRunning = false; end
        end
        function rtn = setControlMode(self, m)
            arguments
                self (1,1) SessionMockController
                m (1,1) {mustBeInteger, mustBeGreaterThanOrEqual(m, 0), mustBeLessThanOrEqual(m, 7)}
            end
            rtn = self.rec('setControlMode', {m});
        end
        function rtn = setPatternID(self, id)
            arguments
                self (1,1) SessionMockController
                id (1,1) {mustBeInteger, mustBeGreaterThanOrEqual(id, 0), mustBeLessThanOrEqual(id, 65535)}
            end
            rtn = self.rec('setPatternID', {id});
        end
        function setPositionX(self, pos)
            arguments
                self (1,1) SessionMockController
                pos (1,1) {mustBeInteger, mustBeGreaterThanOrEqual(pos, 0), mustBeLessThanOrEqual(pos, 65535)}
            end
            self.rec('setPositionX', {pos});
        end
        function rtn = setPatternFunctionID(self, id)
            arguments
                self (1,1) SessionMockController
                id (1,1) {mustBeInteger, mustBeGreaterThanOrEqual(id, 0), mustBeLessThanOrEqual(id, 65535)}
            end
            rtn = self.rec('setPatternFunctionID', {id});
        end
        function rtn = startDisplay(self, deciSeconds, waitForEnd)
            arguments
                self (1,1) SessionMockController
                deciSeconds (1,1) {mustBeInteger, mustBeGreaterThan(deciSeconds, 0), mustBeLessThanOrEqual(deciSeconds, 65535)}
                waitForEnd (1,1) logical = true
            end
            rtn = self.rec('startDisplay', {deciSeconds, waitForEnd});
            if self.sim_speed > 0 && waitForEnd
                pause(self.sim_speed * double(deciSeconds) / 10);
            end
        end
        function rtn = setAOFunctionID(self, ch, id)
            arguments
                self (1,1) SessionMockController
                ch (1,1) {mustBeInteger, mustBeGreaterThanOrEqual(ch, 2), mustBeLessThanOrEqual(ch, 5)}
                id (1,1) {mustBeInteger, mustBeGreaterThanOrEqual(id, 0), mustBeLessThanOrEqual(id, 65535)}
            end
            rtn = self.rec('setAOFunctionID', {ch, id});
        end
        function rtn = setAO(self, ch, voltage)
            arguments
                self (1,1) SessionMockController
                ch (1,1) {mustBeInteger, mustBeGreaterThanOrEqual(ch, 6), mustBeLessThanOrEqual(ch, 7)}
                voltage (1,1) {mustBeNumeric, mustBeFinite, mustBeGreaterThanOrEqual(voltage, -10), mustBeLessThanOrEqual(voltage, 10)}
            end
            rtn = self.rec('setAO', {ch, voltage});
        end

        function names = methodsCalled(self)
            names = cellfun(@(c) c{1}, self.calls, 'UniformOutput', false);
        end

        function t = callTable(self)
            % Call log as a table (method, args as text) for the dry-run report.
            n = numel(self.calls);
            method = cell(n,1); args = cell(n,1);
            for k = 1:n
                method{k} = self.calls{k}{1};
                a = self.calls{k}{2};
                parts = cellfun(@(x) strtrim(formattedDisplayText(x)), a, 'UniformOutput', false);
                args{k} = strjoin(parts, ', ');
            end
            t = table((1:n).', method, args, 'VariableNames', {'call','method','args'});
        end
    end

    methods (Access = private)
        function rtn = rec(self, name, args)
            if ~isfield(self.counts, name), self.counts.(name) = 0; end
            self.counts.(name) = self.counts.(name) + 1;
            self.calls{end+1} = {name, args};
            rtn = true;
            if self.verbose
                fprintf('    [mock] %s(%s)\n', name, strjoin(cellfun(@(x) mat2str(x), args, 'UniformOutput', false), ', '));
            end
            if strcmp(name, self.fail_method) && self.counts.(name) == self.fail_on_call
                if strcmp(self.fail_mode, 'false')
                    rtn = false;
                else
                    error('mock:injected', 'Injected %s failure (call %d)', name, self.counts.(name));
                end
            end
        end
    end
end
