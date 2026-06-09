classdef FilopodiaSegmentationProcess < ImageAnalysisProcess
    % Process 1. Cell body mask + steerable ridge maps (res/theta/nms/scaleMap)
    % used downstream to anchor and trace filopodial shafts.
    % Output base is ImageAnalysisProcess (mixed image-analysis output), unlike
    % FocalAdhesionSegmentationProcess which is mask-only.

    methods (Access = public)
        function obj = FilopodiaSegmentationProcess(owner, varargin)
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
                super_args{2} = FilopodiaSegmentationProcess.getName;
                super_args{3} = @segmentMovieFilopodia;
                if isempty(funParams)
                    funParams = FilopodiaSegmentationProcess.getDefaultParams(owner, outputDir);
                end
                super_args{4} = funParams;
            end
            obj = obj@ImageAnalysisProcess(super_args{:});
        end

        function varargout = loadChannelOutput(obj, iChan, iFrame, varargin)
            outputList = {'bodyMask', 'res', 'theta', 'nms', 'scaleMap'};
            ip = inputParser;
            ip.addRequired('iChan', @(x) isscalar(x) && obj.checkChanNum(x));
            ip.addRequired('iFrame', @(x) isscalar(x) && obj.checkFrameNum(x));
            ip.addParamValue('output', outputList, @(x) all(ismember(x, outputList)));
            ip.addParamValue('useCache', true, @islogical);
            ip.parse(iChan, iFrame, varargin{:});
            output = ip.Results.output;
            if ischar(output), output = {output}; end
            fname = fullfile(obj.outFilePaths_{1, iChan}, ...
                sprintf('filoSeg_frame_%04d.mat', iFrame));
            s = cached.load(fname, '-useCache', ip.Results.useCache, output{:});
            for k = 1:numel(output), varargout{k} = s.(output{k}); end %#ok<AGROW>
        end
    end

    methods (Static)
        function name = getName(), name = 'Filopodia Segmentation'; end
        function h = GUI(), h = @abstractProcessGUI; end   % TODO custom GUI

        function funParams = getDefaultParams(owner, varargin)
            ip = inputParser;
            ip.addRequired('owner', @(x) isa(x, 'MovieData'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            ip.parse(owner, varargin{:});
            outputDir = ip.Results.outputDir;

            funParams.OutputDirectory   = [outputDir filesep 'FilopodiaSegmentation'];
            funParams.ChannelIndex      = 1;            % talin-GFP
            funParams.SteerableOrder    = 4;            % even -> ridge detector
            funParams.SigmaArray        = [1 2 4];      % px; tune to talin PSF
            funParams.BodyThreshold     = 'rosin';      % 'rosin'|'otsu'|scalar
            funParams.BodyMinArea       = 500;          % px
            funParams.BodyClosingRadius = 3;            % px
            funParams.HysteresisHigh    = [];           % []->auto (bkg-based)
            funParams.HysteresisLow     = [];
            funParams.ProcessFrames     = [];           % []->all
            funParams.BatchMode         = false;
        end
    end
end
