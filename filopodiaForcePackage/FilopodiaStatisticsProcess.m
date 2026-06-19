classdef FilopodiaStatisticsProcess < DataProcessingProcess
    % Process 6. Aggregate tracks (P3), windows (P4) and samples (P5) into
    % per-movie statistics: count, length/velocity/fluctuation distributions,
    % lifetime, region-resolved force/talin distributions, and correlations
    % (e.g. tip talin vs tip axial force; protrusion velocity vs tip force).

    methods (Access = public)
        function obj = FilopodiaStatisticsProcess(owner, varargin)
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
                super_args{2} = FilopodiaStatisticsProcess.getName;
                super_args{3} = @computeMovieFilopodiaStats;
                if isempty(funParams)
                    funParams = FilopodiaStatisticsProcess.getDefaultParams(owner, outputDir);
                end
                super_args{4} = funParams;
            end
            obj = obj@DataProcessingProcess(super_args{:});
        end

        function varargout = loadChannelOutput(obj, iChan, varargin)
            outputList = {'stats'};
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


    methods (Access = public)
        function markSuccess(obj)
            % Mark process as successfully completed (for GUI status display).
            % success_ is SetAccess=protected so only subclass methods can set it.
            obj.success_ = true;
            obj.updated_ = true;
        end
    end

    methods (Static)
        function name = getName(), name = 'Filopodia Statistics'; end
        function h = GUI(), h = @filopodiaStatisticsProcessGUI; end

        function funParams = getDefaultParams(owner, varargin)
            ip = inputParser;
            ip.addRequired('owner', @(x) isa(x, 'MovieData'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            ip.parse(owner, varargin{:});
            outputDir = ip.Results.outputDir;

            funParams.OutputDirectory   = [outputDir filesep 'FilopodiaForcePackage' filesep 'FilopodiaStatistics'];
            funParams.ChannelIndex      = 1;
            funParams.ClassProcessIndex = [];     % []->find FilopodiaClassificationProcess (P4)
            funParams.SampleProcessIndex= [];     % []->find FilopodiaSamplingProcess (P5)
            funParams.MinLifetimeForStats = 3;    % frames; ignore very short filopodia
            funParams.ExportCSV         = true;
            funParams.MakeFigures       = true;
            funParams.BatchMode         = false;
        end
    end
end
