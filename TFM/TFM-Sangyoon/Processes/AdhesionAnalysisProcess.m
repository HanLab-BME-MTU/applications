classdef AdhesionAnalysisProcess < Process
    methods (Access = public)
        function obj = AdhesionAnalysisProcess(owner)
            obj = obj@Process(owner, AdhesionAnalysisProcess.getName);
            obj.funName_ = @analyzeAdhesionMaturationProcess;
            obj.funParams_ = AdhesionAnalysisProcess.getDefaultParams(owner);
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
            name = 'AdhesionAnalysis';
        end
        
        function funParams = getDefaultParams(owner)
            % Input check
            ip=inputParser;
            ip.addRequired('owner', @(x) isa(x, 'MovieObject'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            ip.addOptional('iChan',1:numel(owner.channels_),...
               @(x) all(owner.checkChanNum(x)));
            ip.addOptional('showAllTracks',false,@(x)islogical(x)||isempty(x))
            ip.addOptional('plotEachTrack',false,@(x)islogical(x)||isempty(x))
            ip.addOptional('onlyEdge',false,@islogical); % collect NA tracks that ever close to cell edge
            ip.addOptional('outputPath','AdhesionAnalysis',@ischar)
            ip.addOptional('saveAnalysis',true,@islogical)
            ip.addOptional('matchWithFA',true,@islogical) %For cells with only NAs, we turn this off.
            ip.addOptional('minLifetime',5,@isscalar) %For cells with only NAs, we turn this off.
            ip.parse(owner)
            
            % Set default parameters
            funParams.OutputDirectory = [ip.Results.outputDir  filesep 'AdhesionAnalysis'];
            funParams.iChan = ip.Results.iChan;
            funParams.showAllTracks = ip.Results.showAllTracks;
            funParams.plotEachTrack = ip.Results.plotEachTrack;
            funParams.onlyEdge = ip.Results.onlyEdge;
            funParams.saveAnalysis = ip.Results.saveAnalysis;
            funParams.matchWithFA = ip.Results.matchWithFA;
            funParams.minLifetime = ip.Results.minLifetime;
        end
    end
end
