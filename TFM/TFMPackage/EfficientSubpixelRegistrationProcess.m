classdef  EfficientSubpixelRegistrationProcess < StageDriftCorrectionProcess & NonSingularProcess
    % Concrete class for a stage drift correction process
    %
    % Andrew R. Jamieson Feb. 2017
    
    methods
        function obj = EfficientSubpixelRegistrationProcess(owner,varargin)
            
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
                
                % Constructor for EfficientSubpixelRegistrationProcess

                super_args{1} = owner;
                super_args{2} = EfficientSubpixelRegistrationProcess.getName;
                super_args{3} = @efficientSubPixelRegistration;
                if isempty(funParams)
                    funParams = EfficientSubpixelRegistrationProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            obj = obj@StageDriftCorrectionProcess(super_args{:});
        end
        
        
        function output = getDrawableOutput(obj, varargin)
            % Rename default registration output, add bead tracking flow
            output = getDrawableOutput@StageDriftCorrectionProcess(obj);
        end
    end
    methods (Static)
        
        function name = getName()
            name = 'Efficient Subpixel Registration';
        end
        function h = GUI()
            h = @EfficientSubpixelRegistrationProcessGUI;
        end
        
        function funParams = getDefaultParams(owner, varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1 : numel(owner.channels_);
            funParams.OutputDirectory = [outputDir  filesep 'stageDriftCorrection_Efficientsubpixel'];
            funParams.referenceFramePath = '';
            funParams.referenceFrameNum = 1;
            funParams.iBeadChannel = 1;
            funParams.usfac = 20;
        end
    end
end