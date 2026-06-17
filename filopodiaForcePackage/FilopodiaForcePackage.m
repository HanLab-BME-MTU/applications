classdef FilopodiaForcePackage < Package
    % FilopodiaForcePackage  Filopodia mechanics from talin-GFP + bead TFM.
    %
    % Process chain:
    %   1 FilopodiaSegmentationProcess  body mask + steerable ridge maps
    %   2 FilopodiaDetectionProcess     tip/base puncta + anchored shaft trace
    %   3 FilopodiaTrackingProcess      link over time, dL/dt, fluctuation freq
    %   4 FilopodiaClassificationProcess tip/base/shaft labeling + filopodium assembly
    %   5 FilopodiaSamplingProcess      force (ForceFieldCalculationProcess) + talin
    %   6 FilopodiaStatisticsProcess    population statistics & correlations
    %
    % Traction is read from an existing TFMPackage on the same MovieData.
    % Sangyoon J. Han / 2026

    methods
        function obj = FilopodiaForcePackage(owner, varargin)
            if nargin == 0
                super_args = {};
            else
                ip = inputParser;
                ip.addRequired('owner', @(x) isa(x, 'MovieData'));
                ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
                ip.parse(owner, varargin{:});
                outputDir = ip.Results.outputDir;

                super_args{1} = owner;
                super_args{2} = [outputDir filesep 'FilopodiaForcePackage'];
            end
            obj = obj@Package(super_args{:});
        end

        function [status, processExceptions] = sanityCheck(obj, varargin)
            % Movie metadata needed for length/velocity in physical units.
            missing = @(x) sprintf(['Missing %s. Please edit the movie ' ...
                'and fill in the %s.'], x, x);
            assert(~isempty(obj.owner_.pixelSize_),   missing('pixel size'));
            assert(~isempty(obj.owner_.timeInterval_),missing('time interval'));

            [status, processExceptions] = sanityCheck@Package(obj, varargin{:});
        end
    end

    methods (Static)
        function name = getName()
            name = 'Filopodia Force Analysis';
        end

        function varargout = GUI(varargin)
            varargout{1} = packageGUI('FilopodiaForcePackage', varargin{:});
        end

        function classes = getProcessClassNames(index)
            procContrs = {
                'FilopodiaSegmentationProcess', ...
                'FilopodiaDetectionProcess', ...
                'FilopodiaTrackingProcess', ...
                'FilopodiaClassificationProcess', ...
                'FilopodiaSamplingProcess', ...
                'FilopodiaStatisticsProcess'};
            if nargin == 0, index = 1:numel(procContrs); end
            classes = procContrs(index);
        end

        function m = getDependencyMatrix(i, j)
            %    1 2 3 4 5 6  {Processes}
            m = [0 0 0 0 0 0;   % 1 Segmentation
                 1 0 0 0 0 0;   % 2 Detection
                 0 1 0 0 0 0;   % 3 Tracking
                 0 0 1 0 0 0;   % 4 Classification
                 0 0 0 1 0 0;   % 5 Sampling (TFM force is cross-package)
                 0 0 1 1 1 0];  % 6 Statistics
            if nargin < 2, j = 1:size(m, 2); end
            if nargin < 1, i = 1:size(m, 1); end
            m = m(i, j);
        end

        function procConstr = getDefaultProcessConstructors(index)
            procContrs = {
                @(x,y) FilopodiaSegmentationProcess(x, y, FilopodiaForcePackage.getDefaultSegParams(x, y)), ...
                @(x,y) FilopodiaDetectionProcess(x, y, FilopodiaForcePackage.getDefaultDetectionParams(x, y)), ...
                @(x,y) FilopodiaTrackingProcess(x, y, FilopodiaForcePackage.getDefaultTrackingParams(x, y)), ...
                @(x,y) FilopodiaClassificationProcess(x, y, FilopodiaForcePackage.getDefaultClassificationParams(x, y)), ...
                @(x,y) FilopodiaSamplingProcess(x, y, FilopodiaForcePackage.getDefaultSamplingParams(x, y)), ...
                @(x,y) FilopodiaStatisticsProcess(x, y, FilopodiaForcePackage.getDefaultStatsParams(x, y))};
            if nargin == 0, index = 1:numel(procContrs); end
            procConstr = procContrs(index);
        end

        function movieClass = getMovieClass()
            movieClass = 'MovieData';
        end

        % ---- per-process default parameter injectors (mirror FA package) ----
        function funParams = getDefaultSegParams(owner, outputDir)
            funParams = FilopodiaSegmentationProcess.getDefaultParams(owner, outputDir);
            funParams.ChannelIndex = 1;             % talin-GFP
        end
        function funParams = getDefaultDetectionParams(owner, outputDir)
            funParams = FilopodiaDetectionProcess.getDefaultParams(owner, outputDir);
            funParams.ChannelIndex = 1;
            hasPSF = arrayfun(@(x) ~isempty(x.psfSigma_), owner.channels_);
            if any(hasPSF)
                funParams.PSFsigma = owner.channels_(find(hasPSF, 1)).psfSigma_;
            end
        end
        function funParams = getDefaultTrackingParams(owner, outputDir)
            funParams = FilopodiaTrackingProcess.getDefaultParams(owner, outputDir);
        end
        function funParams = getDefaultClassificationParams(owner, outputDir)
            funParams = FilopodiaClassificationProcess.getDefaultParams(owner, outputDir);
        end
        function funParams = getDefaultSamplingParams(owner, outputDir)
            funParams = FilopodiaSamplingProcess.getDefaultParams(owner, outputDir);
            funParams.IntensityChannel = 1;         % talin-GFP
        end
        function funParams = getDefaultStatsParams(owner, outputDir)
            funParams = FilopodiaStatisticsProcess.getDefaultParams(owner, outputDir);
        end
    end
end
