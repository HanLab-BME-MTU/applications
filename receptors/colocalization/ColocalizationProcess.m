classdef ColocalizationProcess < ImageAnalysisProcess
    
        methods (Access = public)
        function obj = ColocalizationProcess(owner,varargin)
            
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
                super_args{2} = ColocalizationProcess.getName;
                super_args{3} = @colocalizationWrapper;
                if isempty(funParams)
                    funParams = ColocalizationProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
                if(nargin > 4)
                    super_args{5:nargin} = varargin{5:nargin};
                end
            end
            
            obj = obj@ImageAnalysisProcess(super_args{:});
        end
        
    end
    methods (Static)
        function name = getName()
            name = 'Colocalization';
        end

        function methods = getMethods(varargin)
            colocalizationMethods(1).name = 'Point2Point';
            colocalizationMethods(1).func = @colocalMeasurePt2Pt;
            colocalizationMethods(2).name = 'Point2Continuum';
            colocalizationMethods(2).func = @colocalMeasurePt2Cnt;
            colocalizationMethods(3).name = 'Continuum2Continuum';
            colocalizationMethods(3).func = @colocalMeasureCnt2Cnt;
            
            ip=inputParser;
            ip.addOptional('index',1:length(colocalizationMethods),@isvector);
            ip.parse(varargin{:});
            index = ip.Results.index;
            methods=colocalizationMethods(index);
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelRef = 1;
            funParams.ChannelObs = 2;
            funParams.ChannelMask = 2;
            funParams.SearchRadius = 3;
            funParams.RandomRuns =1;
            funParams.OutputDirectory = [outputDir  filesep 'colocalInfo.mat'];
            funParams.ProcessIndex = [];%Default is to use raw images
            funParams.MethodIndx = 2;
        end
    end
end