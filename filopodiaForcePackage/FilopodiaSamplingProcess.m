classdef FilopodiaSamplingProcess < DataProcessingProcess
    % Process 5. Per window per frame: sample reconstructed traction
    % (magnitude + axial/lateral via projection on the local tangent) and
    % talin-GFP intensity. Traction is read from the TFMPackage's
    % ForceFieldCalculationProcess output 'tMapUnshifted' (lccb convention).

    methods (Access = public)
        function obj = FilopodiaSamplingProcess(owner, varargin)
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
                super_args{2} = FilopodiaSamplingProcess.getName;
                super_args{3} = @sampleMovieFilopodia;
                if isempty(funParams)
                    funParams = FilopodiaSamplingProcess.getDefaultParams(owner, outputDir);
                end
                super_args{4} = funParams;
            end
            obj = obj@DataProcessingProcess(super_args{:});
        end

        function varargout = loadChannelOutput(obj, iChan, varargin)
            outputList = {'filoSamples'};
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
        function name = getName(), name = 'Filopodia Force / Intensity Sampling'; end
        function h = GUI(), h = @filopodiaSamplingProcessGUI; end

        function funParams = getDefaultParams(owner, varargin)
            ip = inputParser;
            ip.addRequired('owner', @(x) isa(x, 'MovieData'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            ip.parse(owner, varargin{:});
            outputDir = ip.Results.outputDir;

            funParams.OutputDirectory   = [outputDir filesep 'FilopodiaForcePackage' filesep 'FilopodiaSampling'];
            funParams.ChannelIndex      = 1;
            funParams.ClassProcessIndex = [];     % []->find FilopodiaClassificationProcess (P4)
            % cross-package traction source (TFMPackage)
            funParams.ForcePackageName  = 'TFMPackage';
            funParams.ForceProcessName  = 'ForceFieldCalculationProcess';
            % sampling
            funParams.ShaftSampleStep   = 3;       % px; arc-length step along shaft
            funParams.SampleRadius      = 1;       % px; local averaging radius for talin
            funParams.SampleStat        = 'mean';  % 'mean'|'median'|'max'
            funParams.BatchMode         = false;
        end
    end
end
