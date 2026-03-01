classdef OtherChannelSamplingProcess < DataProcessingProcess
    % OtherChannelSamplingProcess
    %
    % Samples intensity statistics from a specified channel inside a binary mask
    % for each frame. Calls sampleMovieOtherChannel().
    %
    % 2026 - Sangyoon Han lab

    methods
        function obj = OtherChannelSamplingProcess(owner, varargin)

            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner', @(x) isa(x,'MovieData'));
                ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
                ip.addOptional('funParams', [], @isstruct);
                ip.parse(owner, varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;

                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = OtherChannelSamplingProcess.getName;
                super_args{3} = @sampleMovieOtherChannel;

                if isempty(funParams)
                    funParams = OtherChannelSamplingProcess.getDefaultParams(owner, outputDir);
                end
                super_args{4} = funParams;
            end

            obj = obj@DataProcessingProcess(super_args{:});
        end

        function status = checkChannelOutput(obj, varargin)
            % Single output .mat file
            status = false;
            try
                status = logical(~isempty(obj.outFilePaths_) && exist(obj.outFilePaths_{1}, 'file'));
            catch
                status = false;
            end
        end
    end

    methods (Static)
        function name = getName()
            name = 'OtherChannelSamplingProcess';
        end
        
        function h = GUI()
            h= @otherChannelSamplingProcessGUI;
        end

        function funParams = getDefaultParams(owner, varargin)
            % Input check
            ip = inputParser;
            ip.addRequired('owner', @(x) isa(x,'MovieData'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            ip.parse(owner, varargin{:});
            outputDir = ip.Results.outputDir;

            % Defaults
            funParams.OutputDirectory = [outputDir filesep 'otherChannelSampling'];

            % Target channel to sample
            funParams.ChannelIndex = 1;

            % Mask source (process class name and channel index)
            funParams.MaskProcessName = '';  % if empty, auto-detect in sampleMovieOtherChannel
            funParams.MaskChannelIndex = 1;

            % Optional dependencies / options
            funParams.UseStageDriftCorrection = false;

            % Save preview images with mask outline
            funParams.SavePerFrameTifPreview = false;

            % Per-object stats using connected components on mask
            funParams.UseLabeling = false;
            funParams.MinAreaPix = 50;

            % dF/F0
            funParams.ComputeDFF0 = false;
            funParams.BaselineFrames = 1;
        end
    end
end
