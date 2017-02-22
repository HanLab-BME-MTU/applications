classdef AdhesionAnalysisProcess < ImageAnalysisProcess
    
    methods (Access = public)
    
        function obj = AdhesionAnalysisProcess(owner, varargin)
    
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
                super_args{2} = AdhesionAnalysisProcess.getName;
                super_args{3} = @analyzeAdhesionMaturation;
                if isempty(funParams)
                    funParams = AdhesionAnalysisProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            obj = obj@ImageAnalysisProcess(super_args{:});
            
            
        end          
            
        function output = loadChannelOutput(obj, iChan, varargin)
            outputList = {};
            nOutput = length(outputList);

            ip.addRequired('iChan',@(x) obj.checkChanNum(x));
            ip.addOptional('iOutput',1,@(x) ismember(x,1:nOutput));
            ip.addParamValue('output','',@(x) all(ismember(x,outputList)));
            ip.addParamValue('useCache',false,@islogical);
            ip.parse(iChan,varargin{:})
    
            s = cached.load(obj.outFilePaths_{iChan},'-useCache',ip.Results.useCache);

            output = s.Imean;          
        end
    end
    methods (Static)
        function name = getName()
            name = 'Focal Adhesion Analysis';
        end

        function h = GUI()
            h = @AdhesionAnalysisProcessGUI;
        end
        
        function funParams = getDefaultParams(owner, varargin)

            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
%             ip.addOptional('iChan', 1:numel(owner.channels_),...
%                @(x) all(owner.checkChanNum(x)));
            ip.parse(owner, varargin{:});
            
            % Set default parameters
            funParams.OutputDirectory = [ip.Results.outputDir filesep 'AdhesionAnalysis'];

            funParams.showAllTracks = false;
            funParams.plotEachTrack = false;
            funParams.onlyEdge = false; 
            funParams.saveAnalysis = true;
            funParams.matchWithFA = true; 
            funParams.minLifetime = 5; 
            funParams.reTrack = true;
            funParams.skipOnlyReading = false;
            funParams.getEdgeRelatedFeatures = true;

            % Make option for using segmentations?
            % --> see field  "existSegmentation= false;""
            
        end
    end
end
