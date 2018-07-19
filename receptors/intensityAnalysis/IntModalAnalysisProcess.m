classdef IntModalAnalysisProcess < DataProcessingProcess
    %Class definition for intensity modal analysis
        
    methods(Access = public)   
        
        function obj = IntModalAnalysisProcess(owner,varargin)
            
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
                super_args{2} = IntModalAnalysisProcess.getName;
                super_args{3} = @analyzeIntensityModesOO;                               
                if isempty(funParams)                                       
                    funParams = IntModalAnalysisProcess.getDefaultParams(owner,outputDir);                               
                end
                super_args{4} = funParams;                    
            end
            
            obj = obj@DataProcessingProcess(super_args{:});
        end        
        
    end
    
    methods(Static)
        
        function name =getName()
            name = 'Intensity Modal Analysis';
        end

        function funParams = getDefaultParams(owner,varargin)
            
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:});
            outputDir=ip.Results.outputDir;
            
            % Define default process parameters
            %             funParams=IntModalAnalysisProcess.getDefaultParams(owner,outputDir);
            funParams.OutputDirectory = [outputDir  filesep 'IntModalAnalysis'];
            funParams.movieName = 'movie1';
            funParams.startFrame = 1;
            funParams.endFrame = [];
            funParams.alpha = 0.05;
            funParams.variableMean = 0;
            funParams.variableStd = 0;
            funParams.numModeMinMax = [1 9];
            funParams.plotResults = 1;
            funParams.logData = 0;
            funParams.modeParamIn = [];
            funParams.ampOrInt = 1;
            funParams.saveFullFileName = [];
            funParams.ratioTol = 0;
            
        end
        
    end
    
end
