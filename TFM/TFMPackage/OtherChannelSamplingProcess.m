
classdef OtherChannelSamplingProcess < DataProcessingProcess
    % OtherChannelSamplingProcess
    % Samples intensity statistics from an arbitrary MovieData channel
    % (e.g., Fluo-4 Ca2+ channel) inside an existing binary mask per frame.
    %
    % Outputs:
    %   - otherChannelSampling.mat (per-frame summary + optional per-object stats)
    %
    % Sangyoon Han lab - added 2026

    methods
        function obj = OtherChannelSamplingProcess(owner, varargin)
            % obj = OtherChannelSamplingProcess(owner, outputDir, funParams)
            ip = inputParser;
            ip.addRequired('owner', @(x) isa(x,'MovieData'));
            ip.addOptional('outputDir', '', @(x) ischar(x) || isstring(x));
            ip.addOptional('funParams', [], @(x) isempty(x) || isstruct(x));
            ip.parse(owner, varargin{:});

            outputDir = char(ip.Results.outputDir);
            funParamsIn = ip.Results.funParams;

            if isempty(outputDir)
                outputDir = fullfile(owner.outputDirectory_, 'otherChannelSampling');
            end

            funParams = OtherChannelSamplingProcess.getDefaultParams(owner, outputDir);
            if ~isempty(funParamsIn)
                fn = fieldnames(funParamsIn);
                for k = 1:numel(fn)
                    funParams.(fn{k}) = funParamsIn.(fn{k});
                end
            end

            obj@DataProcessingProcess(owner, OtherChannelSamplingProcess.getName(), ...
                @sampleMovieOtherChannel, funParams);
        end
    end

    methods (Static)
        function name = getName()
            name = 'Other Channel Sampling';
        end

        function funParams = getDefaultParams(owner, outputDir) %#ok<INUSD>
            funParams.OutputDirectory   = outputDir;

            % Channel to sample (e.g., Ca2+ channel)
            funParams.ChannelIndex      = 1;

            % Mask source:
            %  - if empty, auto-detect MaskIntersectionProcess, MaskRefinementProcess, or ThresholdProcess
            funParams.MaskProcessName   = '';   % e.g., 'MaskRefinementProcess'
            funParams.MaskChannelIndex  = 1;    % which channel the mask corresponds to
            funParams.UseStageDriftCorrection = true; % warp masks using SDC transforms if available

            % Per-object stats (connected components per frame)
            funParams.UseLabeling       = true;
            funParams.MinAreaPix        = 50;

            % dF/F0
            funParams.ComputeDFF0       = true;
            funParams.BaselineFrames    = 1:10;  % used for F0

            % Saving
            funParams.SavePerFrameTifPreview = false; % optional debugging
        end
    end
end
