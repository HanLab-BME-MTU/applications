classdef FocalAdhesionSegmentationProcess < SegmentationProcess
    %A process for segmenting focal adhesions
    %Hunter Elliott
    %3/2013
    
    methods (Access = public)
        function obj = FocalAdhesionSegmentationProcess(owner,varargin)
            
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
                super_args{2} = FocalAdhesionSegmentationProcess.getName;
                super_args{3} = @segmentFocalAdhesions;
                if isempty(funParams)
                    funParams = FocalAdhesionSegmentationProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            obj = obj@SegmentationProcess(super_args{:});
        end
        
    end
    methods (Static)
        function name = getName()
            name = 'Focal Adhesion Segmentation';
        end
        function h = GUI()
            h= @focalAdhesionSegmentationProcesGUI;
        end
%         function methods = getMethods(varargin)
%             thresholdingMethods(1).name = 'MinMax';
%             thresholdingMethods(1).func = @thresholdFluorescenceImage;
%             thresholdingMethods(2).name = 'Otsu';
%             thresholdingMethods(2).func = @thresholdOtsu;
%             thresholdingMethods(3).name = 'Rosin';
%             thresholdingMethods(3).func = @thresholdRosin;
%             thresholdingMethods(4).name = 'Gradient-based';
%             thresholdingMethods(4).func = @intensityBinnedGradientThreshold;
%             
%             ip=inputParser;
%             ip.addOptional('index',1:length(thresholdingMethods),@isvector);
%             ip.parse(varargin{:});
%             index = ip.Results.index;
%             methods=thresholdingMethods(index);
%         end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1:numel(owner.channels_);
            funParams.OutputDirectory = [outputDir  filesep 'adhesion masks'];
            funParams.SteerableFilterSigma = 250;%Sigma in nm of steerable filter to use in splitting adjacent adhesions
            funParams.OpeningRadiusXY = 0; %Spatial radius in nm of structuring element used in opening.
            funParams.OpeningHeightT = 50; %Temporal "height" in seconds of structuring element used in opening            
            funParams.MinVolTime = 5; %Minimum spatiotemporal "Volume" in micron^2 * seconds of segmented adhesions to retain.
            %funParams.ProcessIndex = [];%Default is to use raw images
            funParams.BatchMode = false;
            %funParams.MethodIndx = 1;
        end
    end
end