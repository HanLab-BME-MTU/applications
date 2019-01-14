classdef TractionForceReadingProcess < DataProcessingProcess
    methods (Access = public)
        function obj = TractionForceReadingProcess(owner,varargin)
%             obj = obj@DataProcessingProcess(owner, TractionForceReadingProcess.getName);
%             obj.funName_ = @readTractionForceFromTracks; % This should be variation from colocalizationAdhesionWithTFM
%             obj.funParams_ = TractionForceReadingProcess.getDefaultParams(owner,varargin{1});
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
                super_args{2} = TractionForceReadingProcess.getName;
                super_args{3} = @readTractionForceFromTracks;
                
                if isempty(funParams)
                    funParams = TractionForceReadingProcess.getDefaultParams(owner,outputDir);
                end
                
                super_args{4} = funParams;
            end
            
            obj = obj@DataProcessingProcess(super_args{:});
        end
        
        function output = loadChannelOutput(obj, iChan, varargin)
            outputList = {'tracksForceMag'};
            nOutput = length(outputList);

            ip =inputParser;
            ip.addRequired('obj');
            ip.addRequired('iChan', @(x) obj.checkChanNum(x));
            ip.addOptional('iOutput',1,@(x) ismember(x,1:nOutput));
            ip.addParamValue('output','',@(x) all(ismember(x,outputList)));
            ip.addParamValue('useCache',false,@islogical);
            ip.addParameter('idSelected', [], @(x) isempty(x) || isnumeric(x));
            ip.parse(obj,iChan,varargin{:})
            idSelected = ip.Results.idSelected;
            output = ip.Results.output;
    
%             s = cached.load(obj.outFilePaths_{iChan},'-useCache',ip.Results.useCache);
            if strcmp(output, 'tracksForceMag')
                if isempty(idSelected)
                    iAdhProc = obj.owner_.getProcessIndex('AdhesionAnalysisProcess');
                    adhAnalProc = obj.owner_.getProcess(iAdhProc);
                    s = load(adhAnalProc.outFilePaths_{1,iChan},'metaTrackData');
                    metaTrackData = s.metaTrackData;
                    loadingSequence=metaTrackData.numTracks:-1:1;
                else
                    loadingSequence=idSelected;
                end

                jj=0;
                for ii=loadingSequence
                    if ~isempty(idSelected)
                        jj=jj+1;
                        progressText((jj)/numel(loadingSequence),'Loading tracksNA') % Update text
                    else
                        jj=ii;
                        progressText((numel(loadingSequence)-ii)/numel(loadingSequence),'Loading tracksNA') % Update text
                    end
                    try
                        curTrackObj = load(trackIndPath(ii),'curTrack');
                    catch
                        continue
                    end
                    if ii~=loadingSequence(1)
                        potentialExtraFields = setdiff(fieldnames(curTrackObj.curTrack),fieldnames(tracksNA));
                        if ~isempty(potentialExtraFields)
                            for pp=1:numel(potentialExtraFields)
                                tracksNA(end).(potentialExtraFields{pp})=[];
                            end
                        end
                    end
                    tracksNA(jj,1) = curTrackObj.curTrack;
                end
                % Might need to filter out failed tracks
                indEmptyTracks = arrayfun(@(x) isempty(x.xCoord),tracksNA);
                tracksNA = tracksNA(~indEmptyTracks);
                
                varargout{iout} = tracksNA;                                 
            end
%             output = s.Imean;          
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
