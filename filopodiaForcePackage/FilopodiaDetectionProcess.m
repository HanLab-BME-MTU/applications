classdef FilopodiaDetectionProcess < DetectionProcess
    % Process 2. pointSourceDetection on talin -> bright puncta; classify
    % tip vs base/FA by body mask; pair tip<->base by LAP; trace the dim
    % shaft as a min-cost path on the steerable response between the two
    % anchored endpoints. Saves a tracker-compatible movieInfo plus filoInfo.

    methods (Access = public)
        function obj = FilopodiaDetectionProcess(owner, varargin)
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

                super_args{1} = owner;
                super_args{2} = FilopodiaDetectionProcess.getName;
                super_args{3} = @detectMovieFilopodia;
                if isempty(funParams)
                    funParams = FilopodiaDetectionProcess.getDefaultParams(owner, outputDir);
                end
                super_args{4} = funParams;
            end
            obj = obj@DetectionProcess(super_args{:});
        end

        function varargout = loadChannelOutput(obj, iChan, varargin)
            outputList = {'movieInfo', 'filoInfo'};
            ip = inputParser;
            ip.addRequired('iChan', @(x) isscalar(x) && obj.checkChanNum(x));
            ip.addOptional('iFrame', 1:obj.owner_.nFrames_, @(x) all(obj.checkFrameNum(x)));
            ip.addParamValue('output', outputList, @(x) all(ismember(x, outputList)));
            ip.addParamValue('useCache', true, @islogical);
            ip.parse(iChan, varargin{:});
            iFrame = ip.Results.iFrame; output = ip.Results.output;
            if ischar(output), output = {output}; end
            s = cached.load(obj.outFilePaths_{1, iChan}, '-useCache', ip.Results.useCache, output{:});
            for k = 1:numel(output)
                val = s.(output{k});
                if numel(iFrame) == 1 && iscell(val), val = val{iFrame};
                elseif numel(iFrame) == 1 && isstruct(val) && numel(val) >= iFrame
                    val = val(iFrame);
                end
                varargout{k} = val; %#ok<AGROW>
            end
        end
    end

    methods (Static)
        function name = getName(), name = 'Filopodia Detection (tip/base/shaft)'; end
        function h = GUI(), h = @abstractProcessGUI; end   % TODO custom GUI

        function funParams = getDefaultParams(owner, varargin)
            ip = inputParser;
            ip.addRequired('owner', @(x) isa(x, 'MovieData'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            ip.parse(owner, varargin{:});
            outputDir = ip.Results.outputDir;

            funParams.OutputDirectory = [outputDir filesep 'FilopodiaDetection'];
            funParams.ChannelIndex    = 1;          % talin-GFP
            funParams.SegProcessIndex = [];         % []->find FilopodiaSegmentationProcess
            % pointSourceDetection
            funParams.PSFsigma        = 1.5;        % px; from channel psfSigma_
            funParams.Alpha           = 0.05;
            funParams.ConfRadius      = [];         % default 2*sigma
            funParams.WindowSize      = [];         % default 4*sigma
            % tip vs base classification
            funParams.TipMaxDistFromBody = 50;      % px outside body -> tip candidate
            funParams.BaseSearchBand     = 5;       % px band around body edge -> base
            % shaft tracing
            funParams.MaxTipBaseDist  = 100;        % px, max plausible length
            funParams.OrientTolerance = 30;         % deg
            funParams.OrientLambda    = 2;          % orientation penalty weight
            funParams.PathCostFloor   = 1e-3;
            funParams.MinFiloLength   = 5;          % px
            funParams.ProcessFrames   = [];
            funParams.BatchMode       = false;
        end
    end
end
