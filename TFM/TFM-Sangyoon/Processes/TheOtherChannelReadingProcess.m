classdef TheOtherChannelReadingProcess < DataProcessingProcess
    methods (Access = public)
        function obj = TheOtherChannelReadingProcess(owner,varargin)
%             obj = obj@DataProcessingProcess(owner, TheOtherChannelReadingProcess.getName);
%             obj.funName_ = @readTheOtherChannelFromTracks;
%             obj.funParams_ = TheOtherChannelReadingProcess.getDefaultParams(owner,varargin{1});
            
            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieData'));
                ip.addOptional('outputDir', owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = TheOtherChannelReadingProcess.getName;
                super_args{3} = @readTheOtherChannelFromTracks;
                
                if isempty(funParams)
                    funParams = TheOtherChannelReadingProcess.getDefaultParams(owner,outputDir);
                end
                
                super_args{4} = funParams;
            end
            
            obj = obj@DataProcessingProcess(super_args{:});
        end
        
        function output = loadChannelOutput(obj, iChan, varargin)
            ip =inputParser;
            ip.addRequired('obj');
            ip.addRequired('iChan',@(x) obj.checkChanNum(x));
            ip.addParamValue('useCache',false,@islogical);
            ip.parse(obj,iChan,varargin{:})
    
            % relocate metaTrackData.trackFolderPath with current
            % directory
            outputFile = obj.outFilePaths_;
            [prevProcessPath,trackIndividualName] = fileparts(obj.outFilePaths_{iChan});
            currentProcessPath = [obj.owner_.outputDirectory_ filesep 'FocalAdhesionPackage' filesep 'TheOtherChannelReading'];
            if ~strcmp(prevProcessPath,currentProcessPath)
                outputFile{iChan} = [currentProcessPath filesep trackIndividualName '.mat'];
                obj.setOutFilePaths(outputFile);
            end

            output = cached.load(obj.outFilePaths_{iChan},'-useCache',ip.Results.useCache);
%             output = s.Imean;          
        end
    end
    methods (Static)
        function name = getName()
            name = 'The Other Channel Reading';
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
            ip.addOptional('iChanSlave',setdiff(1:numel(owner.channels_),pAnal.ChannelIndex),...
               @(x) all(owner.checkChanNum(x)));
            ip.addParameter('doFAregistration', true, @islogical);
            ip.parse(owner,varargin{:})
            
            % Set default parameters
            funParams.OutputDirectory = [ip.Results.outputDir filesep 'TheOtherChannelReading'];
            funParams.ChannelIndex = ip.Results.ChannelIndex;
            funParams.iChanSlave = ip.Results.iChanSlave;
            funParams.doFAregistration = ip.Results.doFAregistration;
        end
        
        function h = GUI()
            h = @theOtherChannelReadingProcessGUI;
        end
    end
end
