classdef FilopodiaTrackingProcess < DataProcessingProcess
    % Process 3. Link tip (+base) detections over time by LAP (lap.m from
    % uTrack), gap-close, then derive L(t)=geodesic base->tip, signed dL/dt
    % (protrusion +, retraction -), motility state and fluctuation frequency.
    % To reuse the full uTrack engine instead, swap the wrapper for
    % trackCloseGapsKalmanSparse on the movieInfo from process 2.

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

                super_args{1} = owner;
                super_args{2} = FilopodiaTrackingProcess.getName;
                super_args{3} = @trackMovieFilopodia;
                if isempty(funParams)
                    funParams = FilopodiaTrackingProcess.getDefaultParams(owner, outputDir);
                end
                super_args{4} = funParams;
            end
            obj = obj@DataProcessingProcess(super_args{:});
        end

        function varargout = loadChannelOutput(obj, iChan, varargin)
            outputList = {'filoTracks'};
            ip = inputParser;
            ip.addRequired('iChan', @(x) isscalar(x) && obj.checkChanNum(x));
            ip.addParamValue('output', outputList, @(x) all(ismember(x, outputList)));
            ip.addParamValue('useCache', true, @islogical);
            ip.parse(iChan, varargin{:});
            output = ip.Results.output; if ischar(output), output = {output}; end
            s = cached.load(obj.outFilePaths_{1, iChan}, '-useCache', ip.Results.useCache, output{:});
            for k = 1:numel(output), varargout{k} = s.(output{k}); end %#ok<AGROW>
        end
    end

    methods (Static)
        function name = getName(), name = 'Filopodia Tracking'; end
        function h = GUI(), h = @abstractProcessGUI; end

        function funParams = getDefaultParams(owner, varargin)
            ip = inputParser;
            ip.addRequired('owner', @(x) isa(x, 'MovieData'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            ip.parse(owner, varargin{:});
            outputDir = ip.Results.outputDir;

            funParams.OutputDirectory = [outputDir filesep 'FilopodiaForcePackage' filesep 'FilopodiaTracking'];
            funParams.ChannelIndex    = 1;
            funParams.DetProcessIndex = [];      % []->find FilopodiaDetectionProcess
            funParams.MaxLinkDist     = 10;      % px / frame (tip)
            funParams.LinkUseBase     = true;
            funParams.MaxGapFrames    = 2;
            funParams.MinTrackLength  = 3;       % frames
            funParams.VelSmoothWin    = 3;       % frames
            funParams.PauseThreshVel  = [];      % length/frame; []->auto
            funParams.FreqMethod      = 'psd';   % 'psd'|'autocorr'|'cyclecount'
            funParams.DetrendL        = true;
            funParams.BatchMode       = false;
        end
    end
end
