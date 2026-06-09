classdef FilopodiaWindowingProcess < DataProcessingProcess
    % Process 4. Lay temporally consistent windows along each tracked
    % filopodium centerline, parameterized by arclength (base s=0 -> tip s=L).
    % Stores region label, arclength, local tangent and two footprints:
    % a narrow one for talin intensity and a wider one for TFM-resolution force.

    methods (Access = public)
        function obj = FilopodiaWindowingProcess(owner, varargin)
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
                super_args{2} = FilopodiaWindowingProcess.getName;
                super_args{3} = @windowMovieFilopodia;
                if isempty(funParams)
                    funParams = FilopodiaWindowingProcess.getDefaultParams(owner, outputDir);
                end
                super_args{4} = funParams;
            end
            obj = obj@DataProcessingProcess(super_args{:});
        end

        function varargout = loadChannelOutput(obj, iChan, varargin)
            outputList = {'filoWindows'};
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
        function name = getName(), name = 'Filopodia Windowing'; end
        function h = GUI(), h = @abstractProcessGUI; end

        function funParams = getDefaultParams(owner, varargin)
            ip = inputParser;
            ip.addRequired('owner', @(x) isa(x, 'MovieData'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            ip.parse(owner, varargin{:});
            outputDir = ip.Results.outputDir;

            funParams.OutputDirectory  = [outputDir filesep 'FilopodiaWindowing'];
            funParams.ChannelIndex     = 1;
            funParams.TrackProcessIndex = [];     % []->find FilopodiaTrackingProcess
            funParams.WindowMode       = 'normalized';  % 'normalized'|'fixed'
            funParams.NumWindows       = 10;      % normalized: constant count
            funParams.WindowLength     = 5;       % fixed: px per window
            funParams.LateralWidth     = 3;       % px, intensity footprint
            funParams.ForceLateralWidth = [];     % px, []->auto from TFM reg scale
            funParams.BaseFraction     = 0.15;    % s<=0.15 -> base
            funParams.TipFraction      = 0.15;    % s>=0.85 -> tip
            funParams.TipAnchor        = 'distalEndpoint';
            funParams.BatchMode        = false;
        end
    end
end
