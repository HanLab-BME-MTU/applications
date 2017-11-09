classdef TheOtherChannelReadingProcess < Process
    methods (Access = public)
        function obj = TheOtherChannelReadingProcess(owner)
            obj = obj@Process(owner, TheOtherChannelReadingProcess.getName);
            obj.funName_ = @readTheOtherChannelFromTracks;
            obj.funParams_ = TheOtherChannelReadingProcess.getDefaultParams(owner);
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
            name = 'TheOtherChannelReading';
        end
        
        function funParams = getDefaultParams(owner)
            % Input check
            ip=inputParser;
            ip.addRequired('owner', @(x) isa(x, 'MovieObject'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            ip.addOptional('iChanMaster',1,...
               @(x) all(owner.checkChanNum(x)));
            ip.addOptional('iChanSlave',min(2,numel(owner.channels_)),...
               @(x) all(owner.checkChanNum(x)));
            ip.addOptional('outputPath','analysis1',@ischar)
            ip.parse(owner)
            
            % Set default parameters
            funParams.OutputDirectory = [ip.Results.outputDir filesep 'TheOtherChannelReading'];
            funParams.outputPath = ip.Results.outputPath; %This is a specific folder
            funParams.iChanMaster = ip.Results.iChanMaster;
            funParams.iChanSlave = ip.Results.iChanSlave;
            funParams.onlyEdge = ip.Results.onlyEdge;
        end
    end
end
