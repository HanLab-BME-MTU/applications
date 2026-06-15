classdef FilopodiaSegmentationProcess < ImageAnalysisProcess

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
            % Load one frame's segmentation result directly from the output
            % directory stored in funParams_ (robust to outFilePaths_ being empty).
            outputList = {'bodyMask','res','theta','nms','scaleMap','shaftMask'};
            ip = inputParser;
            ip.addRequired('iChan',  @(x) isscalar(x) && obj.checkChanNum(x));
            ip.addRequired('iFrame', @(x) isscalar(x) && obj.checkFrameNum(x));
            ip.addParamValue('output',   outputList, @(x) all(ismember(x, outputList)));
            ip.addParamValue('useCache', true,       @islogical);
            ip.parse(iChan, iFrame, varargin{:});
            output = ip.Results.output;
            if ischar(output), output = {output}; end

            outDir = obj.funParams_.OutputDirectory;
            fname  = fullfile(outDir, sprintf('filoSeg_frame_%04d.mat', iFrame));
            assert(exist(fname,'file')==2, ...
                'P1 frame file not found: %s\nRun FilopodiaSegmentationProcess first.', fname);
            s = cached.load(fname, '-useCache', ip.Results.useCache, output{:});
            for k = 1:numel(output), varargout{k} = s.(output{k}); end %#ok<AGROW>
        end
    end

    methods (Static)
        function name = getName(), name = 'Filopodia Segmentation'; end
        function h = GUI(), h = @abstractProcessGUI; end

        function funParams = getDefaultParams(owner, varargin)
            ip = inputParser;
            ip.addRequired('owner', @(x) isa(x, 'MovieData'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            ip.parse(owner, varargin{:});
            outputDir = ip.Results.outputDir;

            funParams.OutputDirectory   = [outputDir filesep 'FilopodiaForcePackage' filesep 'FilopodiaSegmentation'];
            funParams.ChannelIndex      = 1;
            funParams.SteerableOrder    = 4;
            funParams.SigmaArray        = [1 2 4];
            funParams.BodyThreshold     = 'otsu';
            funParams.GaussianBlurSigma = 2;     % px; blur before body threshold
            funParams.BodyMinArea       = 500;
            funParams.BodyOpenRadius    = 8;     % px; opening removes filopodia roots (despike)
            funParams.BodyClosingRadius = 8;     % px; closing rounds/smooths the body edge
            funParams.HysteresisHigh    = [];
            funParams.HysteresisLow     = [];
            funParams.ProcessFrames     = [];
            funParams.BatchMode         = false;
        end
    end
end