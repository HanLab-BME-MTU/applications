classdef MultiThreshProcess < ThresholdProcess
    
    methods (Access = public)
        function obj = ThresholdProcess(owner,varargin)
            
            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieData'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = MultiThreshProcess.getName;
                super_args{3} = @thresholdMovie;
                if isempty(funParams)
                    funParams = ThresholdProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            obj = obj@SegmentationProcess(super_args{:});
        end
        
    end
    methods (Static)
        function name = getName()
            name = 'Thresholding';
        end
        function h = GUI()
            h= @thresholdProcessGUI;
        end
        function methods = getMethods(varargin)
            thresholdingMethods(1).name = 'MinMax';
            thresholdingMethods(1).func = @thresholdFluorescenceImage;
            thresholdingMethods(2).name = 'Otsu';
            thresholdingMethods(2).func = @thresholdOtsu;
            thresholdingMethods(3).name = 'Rosin';
            thresholdingMethods(3).func = @thresholdRosin;
            thresholdingMethods(4).name = 'Gradient-based';
            thresholdingMethods(4).func = @intensityBinnedGradientThreshold;
            
            ip=inputParser;
            ip.addOptional('index',1:length(thresholdingMethods),@isvector);
            ip.parse(varargin{:});
            index = ip.Results.index;
            methods=thresholdingMethods(index);
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1:numel(owner.channels_);
            funParams.OutputDirectory = [outputDir  filesep 'masks'];
            funParams.ProcessIndex = [];%Default is to use raw images
            funParams.ThresholdValue = []; % automatic threshold selection
            funParams.MaxJump = 0; %Default is no jump suppression
            funParams.GaussFilterSigma = 0; %Default is no filtering.
            funParams.BatchMode = false;
            funParams.MethodIndx = 1;
        end
    end
end