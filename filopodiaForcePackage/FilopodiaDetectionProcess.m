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

            funParams.OutputDirectory = [outputDir filesep 'FilopodiaForcePackage' filesep 'FilopodiaDetection'];
            funParams.ChannelIndex    = 1;          % talin-GFP
            funParams.SegProcessIndex = [];         % []->find FilopodiaSegmentationProcess
            % pointSourceDetection
            funParams.PSFsigma        = 1.5;        % px; from channel psfSigma_
            funParams.Alpha           = 0.05;
            funParams.ConfRadius      = [];         % default 2*sigma
            funParams.WindowSize      = [];         % default 4*sigma
            % detection mode
            funParams.DetectMode = 'auto';          % 'auto' (multi-frame->all, 1-frame->tip) | 'all' | 'tip'
            funParams.BaseInsideBand = 4;           % px inside body edge still counted as a (base) adhesion ('all' mode)
            % tip vs base classification
            funParams.TipMaxDistFromBody = 130;     % px outside body -> tip/adhesion candidate gate
            funParams.BaseSearchBand     = 5;       % px band around body edge -> base
            funParams.TipRidgeBand       = -1;      % px around shaftMask for ridge gate (<0 = off; cleanup is temporal in P3)
            % shaft tracing
            funParams.MaxTipBaseDist  = 160;        % px, max trace length (~1.2 x TipMaxDistFromBody)
            funParams.OrientTolerance = 30;         % deg
            funParams.OrientLambda    = 3;          % orientation penalty; higher = straighter
            funParams.PathCostFloor   = 1e-3;
            funParams.MinFiloLength   = 5;          % px
            funParams.ShaftAbsorbRadius = 3;        % px; carve width around traced shaft
            funParams.CarveDistalFrac   = 0.8;      % block distal 80% of each shaft
            funParams.ProcessFrames   = [];
            funParams.BatchMode       = false;
        end
    end
end
