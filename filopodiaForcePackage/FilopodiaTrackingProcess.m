classdef FilopodiaTrackingProcess < TrackingProcess
    % Process 3 of the FilopodiaForcePackage.
    %
    % Thin subclass of u-track's TrackingProcess: the full u-track engine
    % (trackCloseGapsKalmanSparse with Brownian + directed-motion cost
    % matrices, gap closing and Kalman filtering) does the linking, and the
    % stock trackingProcessGUI (with its frame-to-frame linking and gap
    % closing setting dialogs) is inherited unchanged. We only (1) seed
    % filopodia-appropriate default parameters and (2) after tracking, convert
    % the u-track tracksFinal into the adhesionTracks struct that Processes
    % 4-6 consume, saving both into the same output file.
    %
    % Input is the FilopodiaDetectionProcess (a DetectionProcess subclass),
    % whose movieInfo is already in u-track format, so no input adapter is
    % needed.

    methods (Access = public)
        function obj = FilopodiaTrackingProcess(owner, varargin)
            if nargin == 0
                super_args = {};
            else
                ip = inputParser;
                ip.addRequired('owner', @(x) isa(x, 'MovieData'));
                ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
                ip.addOptional('funParams', [], @isstruct);
                ip.parse(owner, varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                if isempty(funParams)
                    funParams = FilopodiaTrackingProcess.getDefaultParams(owner, outputDir);
                end
                super_args{1} = owner;
                super_args{2} = outputDir;
                super_args{3} = funParams;
            end
            obj = obj@TrackingProcess(super_args{:});
        end

        function run(obj, varargin)
            % Run the stock u-track tracking, then append adhesionTracks.
            run@TrackingProcess(obj, varargin{:});
            try
                appendAdhesionTracks(obj);
            catch ME
                warning(safeId(ME), 'adhesionTracks conversion failed: %s', ME.message);
            end
        end

        function varargout = loadChannelOutput(obj, iChan, varargin)
            % Extend the parent loader so 'adhesionTracks' is also available
            % to Processes 4-6, while 'tracksFinal' etc. behave as in u-track.
            wantAdh = false;
            for k = 1:numel(varargin)
                if (ischar(varargin{k}) && strcmp(varargin{k},'output')) && k<numel(varargin)
                    o = varargin{k+1};
                    if (ischar(o) && strcmp(o,'adhesionTracks')) || ...
                       (iscell(o) && any(strcmp(o,'adhesionTracks')))
                        wantAdh = true;
                    end
                end
            end
            if wantAdh
                s = load(obj.outFilePaths_{1, iChan}, 'adhesionTracks');
                varargout{1} = s.adhesionTracks;
                return;
            end
            [varargout{1:nargout}] = loadChannelOutput@TrackingProcess(obj, iChan, varargin{:});
        end
    end

    methods (Static)
        function name = getName(), name = 'Filopodia Tracking'; end
        % GUI() inherited from TrackingProcess -> stock trackingProcessGUI

        function funParams = getDefaultParams(owner, varargin)
            ip = inputParser;
            ip.addRequired('owner', @(x) isa(x, 'MovieData'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            ip.parse(owner, varargin{:});
            outputDir = ip.Results.outputDir;

            % start from the full u-track default (gapCloseParam, costMatrices,
            % kalmanFunctions all populated) and adjust for filopodia tips.
            funParams = TrackingProcess.getDefaultParams(owner, outputDir);
            funParams.OutputDirectory = [outputDir filesep 'FilopodiaForcePackage' filesep 'FilopodiaTracking'];

            % point the tracker at our detection process specifically
            funParams.DetProcessIndex = owner.getProcessIndex('FilopodiaDetectionProcess',1,0);

            % filopodia-appropriate gap closing
            funParams.gapCloseParam.timeWindow  = 4;   % max gap (frames)
            funParams.gapCloseParam.mergeSplit  = 0;   % tips don't merge/split
            funParams.gapCloseParam.minTrackLen = 3;   % drop 1-2 frame flickers
            funParams.gapCloseParam.diagnostics = 0;

            % single talin channel by default (not all channels)
            funParams.ChannelIndex = funParams.DetProcessIndex;
            if isempty(funParams.ChannelIndex), funParams.ChannelIndex = 1; end
            funParams.ChannelIndex = 1;
        end
    end
end

% =====================================================================
function appendAdhesionTracks(obj)
% Convert u-track tracksFinal into adhesionTracks and save into the same
% per-channel output file, so Processes 4-6 can read it directly.
MD = obj.getOwner();
pix = MD.pixelSize_; dt = MD.timeInterval_;
chans = obj.funParams_.ChannelIndex;
for i = chans(:)'
    f = obj.outFilePaths_{1, i};
    if isempty(f) || exist(f,'file')~=2, continue; end
    S = load(f, 'tracksFinal');
    adhesionTracks = tracksFinal2adhesionTracks(S.tracksFinal, pix, dt); %#ok<NASGU>
    save(f, 'adhesionTracks', '-append');
end
end

function id = safeId(ME)
id = ME.identifier; if isempty(id), id = 'FilopodiaTrackingProcess:convert'; end
end