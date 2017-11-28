classdef TractionForceReadingProcess < DataProcessingProcess
    methods (Access = public)
        function obj = TractionForceReadingProcess(owner,varargin)
            obj = obj@DataProcessingProcess(owner, TractionForceReadingProcess.getName);
            obj.funName_ = @readTractionForceFromTracks; % This should be variation from colocalizationAdhesionWithTFM
            obj.funParams_ = TractionForceReadingProcess.getDefaultParams(owner,varargin{1});
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
            name = 'Traction Force Reading Process';
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner', @(x) isa(x, 'MovieObject'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            adhAnalProc = owner.getProcess(owner.getProcessIndex('AdhesionAnalysisProcess'));
            pAnal=adhAnalProc.funParams_;
            
            ip.addOptional('ChannelIndex',pAnal.ChannelIndex,...
               @(x) all(owner.checkChanNum(x)));
            ip.parse(owner,varargin{:})
            
            % Set default parameters
            funParams.OutputDirectory = [ip.Results.outputDir filesep 'TractionForceReadingProcess'];
            funParams.ChannelIndex = ip.Results.ChannelIndex;
            funParams.saveTractionField = true;
        end
        
        function h = GUI()
            h = @tractionForceReadingProcessGUI;
        end
    end
end
