classdef AdhesionAnalysisProcess < DataProcessingProcess %& DataProcessingProcess
    
    methods (Access = public)
    
        function obj = AdhesionAnalysisProcess(owner, varargin)
    
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
                super_args{2} = AdhesionAnalysisProcess.getName;
                super_args{3} = @analyzeAdhesionMaturation;
                
                if isempty(funParams)
                    funParams = AdhesionAnalysisProcess.getDefaultParams(owner,outputDir);
                end
                
                super_args{4} = funParams;
            end
            
            obj = obj@DataProcessingProcess(super_args{:});
                        
        end

        function sanityCheck(obj)
            
            sanityCheck@DataProcessingProcess(obj);
            
            % Cell Segmentation Check
            if obj.funParams_.ApplyCellSegMask                
                iProc = obj.funParams_.SegCellMaskProc;
                % Check Mask is available
                if ~isempty(iProc)
                    assert(iProc < length(obj.owner_.processes_), 'Invalid Process # for Cell Mask Process');
                    assert(isa(obj.owner_.getProcess(iProc), 'MaskRefinementProcess'), ['Process: ' num2str(iProc) ' not a MaskRefinementProcess!']);
                    maskProc = obj.owner_.getProcess(iProc);
                    assert(maskProc.checkChannelOutput(obj.funParams_.ChannelIndex), 'Cell Segmentation Mask Output Not found');
                    % iProc = obj.owner_.getProcessIndex('MaskProcess', 'askUser', false, 'nDesired', Inf);
                else
                    iProc = obj.owner_.getProcessIndex('MaskRefinementProcess');   
                    assert(~isempty(iProc), 'MaskRefinementProcess Process not found cannot Apply Cell Mask!');
                    disp('Setting Cell Segmentation Mask Process index');
                    obj.funParams_.SegCellMaskProc = iProc;
                end
            else
                warning('You do not have segmentation process run for a mask. Using entire image field ...')
            end
            
            %% Sanity check for Detecting FAs 
            iProc = obj.funParams_.detectedNAProc;
            if ~isempty(iProc)
                assert(iProc < length(obj.owner_.processes_), ['Invalid Process #' num2str(iProc) 'for FA Detection Process']);
                assert(isa(obj.owner_.getProcess(iProc), 'DetectionProcess'));
            else
                iProc = obj.owner_.getProcessIndex('DetectionProcess');
                assert(~isempty(iProc), 'FA DetectionProcess Process not found!');
                disp('Setting FA DetectionProcess index');
                obj.funParams_.detectedNAProc = iProc;
            end
            
            %% Sanity check for Tracking FAs 
            iProc = obj.funParams_.trackFAProc;
            if ~isempty(iProc)
                assert(iProc < length(obj.owner_.processes_), ['Invalid Process #' num2str(iProc) 'for FA Tracking Process']);
                assert(isa(obj.owner_.getProcess(iProc), 'TrackingProcess'));
            else
                iProc = obj.owner_.getProcessIndex('TrackingProcess');
                assert(~isempty(iProc), 'FA TrackingProcess Process not found!');
                disp('Setting FA TrackingProcess index');
                obj.funParams_.trackFAProc = iProc;
            end
            
            %% Add Sanity check for FAsegmentationProcess Mask
            iProc = obj.funParams_.FAsegProc;
            if ~isempty(iProc)
                assert(iProc < length(obj.owner_.processes_), ['Invalid Process #' num2str(iProc) 'for FA Seg Process']);
                assert(isa(obj.owner_.getProcess(iProc), 'FocalAdhesionSegmentationProcess'));
            else
                iProc = obj.owner_.getProcessIndex('FocalAdhesionSegmentationProcess');
                assert(~isempty(iProc), 'FocalAdhesionSegmentationProcess Process not found!');
                disp('Setting FocalAdhesionSegmentationProcess index');
                obj.funParams_.FAsegProc = iProc;
            end
        end
        
        function varargout = loadChannelOutput(obj, iChan, varargin)
            % Input check
            outputList = {'trackFC','trackNA','trackFA','detBA',...
                          'detectedFA','detFA','detFC','detNA',...
                          'adhboundary_FA', 'adhboundary_FC','tracksNA'};

            ip =inputParser;
            ip.addRequired('obj');
            ip.addRequired('iChan', @(x) obj.checkChanNum(x));
            ip.addOptional('iFrame', [] ,@(x) obj.checkFrameNum(x));
            ip.addParameter('useCache',true,@islogical);
            ip.addParameter('output', outputList{3}, @(x) all(ismember(x,outputList)));
            ip.addParameter('wantFullTrack', false, @islogical);
            ip.addParameter('idSelected', [], @(x) isempty(x) || isnumeric(x));
            ip.parse(obj,iChan,varargin{:})
            output = ip.Results.output;
            varargout = cell(numel(output), 1);
            iFrame = ip.Results.iFrame;
            wantFullTrack = ip.Results.wantFullTrack;
            idSelected = ip.Results.idSelected;
            if ischar(output),output={output}; end
            
            % Data loading
            persistent xCoord yCoord refineFAID stateAll startingFrameExtra endingFrameExtra lastFinishTime 
            if isempty(lastFinishTime)
                lastFinishTime = clock; % assigning current time.. This will be definitely different from obj.finishTime_
            end
            if isempty(xCoord) || isempty(yCoord) || isempty(refineFAID) ...
                    || isempty(stateAll) || isempty(startingFrameExtra) || ...
                    isempty(endingFrameExtra) || ~all(obj.finishTime_==lastFinishTime) ...
                    || strcmp(output,outputList(end))
                try
                    s = load(obj.outFilePaths_{1,iChan},'metaTrackData');
                    metaTrackData = s.metaTrackData;
                    fString = ['%0' num2str(floor(log10(metaTrackData.numTracks))+1) '.f'];
                    numStr = @(trackNum) num2str(trackNum,fString);
                    % relocate metaTrackData.trackFolderPath with current
                    % directory
                    [prevProcessPath,trackIndividualName] = fileparts(metaTrackData.trackFolderPath);
                    currentProcessPath = fileparts(obj.outFilePaths_{1,iChan});
                    if ~strcmp(prevProcessPath,currentProcessPath)
                        metaTrackData.trackFolderPath=[currentProcessPath filesep trackIndividualName];
                        save(obj.outFilePaths_{1,iChan},'metaTrackData');
                    end
                    trackIndPath = @(trackNum) [metaTrackData.trackFolderPath filesep 'track' numStr(trackNum) '.mat'];
                    if isempty(idSelected)
                        loadingSequence=metaTrackData.numTracks:-1:1;
                    else
                        loadingSequence=idSelected;
                    end
                    
                    jj=0;
                    if size(loadingSequence,1)>1
                        loadingSequence=loadingSequence';
                    end
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
                    if strcmp(output,outputList(end))
                        varargout{1} = tracksNA;
                        return
                    end
                catch
                    % Check if the outFilePath has tableTracksNA
                    disp('Checking if the outFilePath has tableTracksNA...')
                    s = load(obj.outFilePaths_{1,iChan},'tracksNA');
                    if isfield(s,'tracksNA')
                        disp('Found the old format. Resaving this with the new format...')
                        % Saving with each track
                        tracksNA = s.tracksNA;
                        trackFolderPath = [obj.funParams_.OutputDirectory filesep 'trackIndividual'];
                        mkdir(trackFolderPath)
                        numTracks = numel(tracksNA);
                        fString = ['%0' num2str(floor(log10(numTracks))+1) '.f'];
                        numStr = @(trackNum) num2str(trackNum,fString);
                        trackIndPath = @(trackNum) [trackFolderPath filesep 'track' numStr(trackNum) '.mat'];

                        for ii=1:numTracks
                            curTrack = tracksNA(ii);
                            save(trackIndPath(ii),'curTrack')
                        end
                        % Saving the metaTrackData
                        metaTrackData.numTracks = numTracks;
                        metaTrackData.trackFolderPath = trackFolderPath;
                        metaTrackData.eachTrackName = 'curTrack';
                        metaTrackData.fString = ['%0' num2str(floor(log10(numTracks))+1) '.f'];
                        metaTrackData.numStr = @(trackNum) num2str(trackNum,fString);
                        metaTrackData.trackIndPath = @(trackNum) [trackFolderPath filesep 'track' numStr(trackNum) '.mat'];
                        save(obj.outFilePaths_{1,iChan},'metaTrackData')
                    end
                end
                if numel(tracksNA)==1
                    s = struct2table(tracksNA,'AsArray',true);
                else
                    s = struct2table(tracksNA);
                end
                xCoord = s.xCoord;
                yCoord = s.yCoord;
                startingFrameExtra = s.startingFrameExtra;
                endingFrameExtra = s.endingFrameExtra;
                refineFAID_cell = s.refineFAID;
                stateAll = s.state;
                % refineFAID_cell is numTracks x
                % refineID_for_everyFramesInvolved. So for each track (each
                % row), I'll make each raw a full frame entries although it
                % is a bit memory intensive
                try %if numel(tracksNA)==1
                    refineFAID = refineFAID_cell;
                catch %else
                    maxFrame = max(cellfun(@length,refineFAID_cell));
                    insuffRows = cellfun(@(x) length(x)<maxFrame,refineFAID_cell);
                    for k=find(insuffRows')
                        refineFAID_cell{k} = [refineFAID_cell{k} ...
                                    NaN(1,maxFrame-length(refineFAID_cell{k}))];
                    end
                    refineFAID = cell2mat(refineFAID_cell);
                end
                lastFinishTime = obj.finishTime_;
            end
            
%             s = cached.load(obj.outFilePaths_{1,iChan}, '-useCache', ip.Results.useCache, 'tableTracksNA');
%             st = cached.load(obj.outFilePaths_{1,iChan}, '-useCache', ip.Results.useCache, 'tracksNA');
            % Note, could do a stack.
            
            %% Check struct vs table loading           
%             if isstruct(s)
%                 s = s.tableTracksNA;
%             else
%                 disp('loaded as table');
%             end
%             s = struct2table(tracksNA);

            nTracks = length(xCoord(:,iFrame));
            number = (1:length(xCoord(:,iFrame)))';
            state = stateAll(:,iFrame); %state = categorical(s.state(:,iFrame));
            iiformat = ['%.' '3' 'd'];
            
            for iout = 1:numel(output)
                switch output{iout}
                    case 'detectedFA'  
                        varargout{1} = t;
                    case 'detBA' 
                        validState = state == 1; %'BA';
                    case {'detNA', 'trackNA'}
                        validState = state == 2; %'NA';
                    case {'detFC', 'trackFC', 'adhboundary_FC'}
                        validState = state == 3; %'FC';
                    case {'detFA', 'trackFA', 'adhboundary_FA'}
                        validState = state == 4; %'FA';
                    case 'staticTracks'
                    case 'tracksNA'

                    otherwise
                        error('Incorrect Output Var type');
                end   
                if strcmp(output{iout}, 'tracksNA')
                
                    varargout{iout} = tracksNA;                                 
                
                elseif ~isempty(strfind(output{iout}, 'det'))
                
                    t = table(xCoord(:,iFrame), yCoord(:,iFrame));
                    varargout{iout} = t{validState,:};                                 
                
                elseif ~isempty(strfind(output{iout},'track'))
                    if ~wantFullTrack
                        vars = {'xCoord', 'yCoord', 'number'};
                        validTracks = validState & startingFrameExtra <= iFrame & endingFrameExtra >= iFrame;                    
                        st = table(xCoord(:,1:iFrame), yCoord(:,1:iFrame), number, ...
                                   'VariableNames', {'xCoord', 'yCoord', 'number'});                    

                        varargout{iout}(nTracks, 1) = struct('xCoord', [], 'yCoord', [], 'number', []);
                        varargout{iout}(validTracks, :) = table2struct(st(validTracks, vars));
                    else
                        varargout{iout} = tracksNA;                                 
                    end
                
                elseif ~isempty(strfind(output{iout},'adhboundary'))                    
                
%                     adhBoundary = cellfun(@(x) x{iFrame}, s{validState, 'adhBoundary'}, 'UniformOutput', false);                         
                    p=obj.funParams_;
                    labelTifPath = [p.OutputDirectory filesep 'labelTifs'];
                    maskAdhesion = imread(strcat(labelTifPath,'/label',num2str(iFrame,iiformat),'.tif'));
                    labelAdhesion = bwlabel(maskAdhesion,4);
                    maxLabel=max(labelAdhesion(:));
                    adhBound = cell(maxLabel,1);
                    for ii=1:maxLabel
                        curAdhBound = bwboundaries(labelAdhesion==ii,4,'noholes');
                        adhBound{ii} = curAdhBound{1}; % strongly assumes each has only one boundary
                    end
                    validAdhState = refineFAID(validState,iFrame); %cellfun(@(x) x(iFrame),refineFAID(validState));
                    
                    varargout{iout} = adhBound(validAdhState); %table2struct(table(adhBoundary, number(validState),'VariableNames',{'adhBoundary', 'number'}));                                 
                else
                    varargout{iout} = [];
                end
            end
            disp(' ')
        end      
    end


    methods (Static)
        function name = getName()
            name = 'Focal Adhesion Analysis';
        end

        function h = GUI()
            h = @focalAdhesionAnalysisProcessGUI;
        end
        
        function output = getDrawableOutput()
            ii = 10;
%             i = ii-1; output(i).name='Before Adhesion Detection'; 
%             output(i).var='detBA';
%             output(i).formatData=[];
%             output(i).type='overlay';
%             output(i).defaultDisplayMethod=@(x) LineDisplay('Marker','.',...
%                 'LineStyle', 'none', 'LineWidth', .6, 'Color', 'g',...
%                 'MarkerSize', 5);            
%             i = ii-2; output(i).name='Nascent Adhesion Detection'; 
%             output(i).var='detNA';
%             output(i).formatData=[];
%             output(i).type='overlay';
%             output(i).defaultDisplayMethod=@(x) LineDisplay('Marker','o',...
%                 'LineStyle','none', 'LineWidth', .8, 'Color', 'r'); 
%             i = ii-3; output(i).name='Focal Contact Detection'; 
%             output(i).var='detFC';
%             output(i).formatData=[];
%             output(i).type='overlay';
%             output(i).defaultDisplayMethod=@(x) LineDisplay('Marker','o',...
%                 'LineStyle','none', 'LineWidth', .8, 'Color', [255/255 153/255 51/255]); 
%             i = ii-4; output(i).name='Focal Adhesion Detection'; 
%             output(i).var='detFA';
%             output(i).formatData=[];
%             output(i).type='overlay';
%             output(i).defaultDisplayMethod=@(x) LineDisplay('Marker','o',...
%                 'LineStyle','none', 'LineWidth', .8, 'Color', 'b'); 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Tracks Display
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            i = ii-5; output(i).name='Nascent Adhesion Tracks'; 
            output(i).var='trackNA';
            output(i).formatData=[];
            output(i).type='overlay';
            output(i).defaultDisplayMethod=@(x) FATracksDisplay('Linewidth', 1, 'Color', 'r'); 
            i = ii-6; output(i).name='Focal Contact Tracks'; 
            output(i).var='trackFC';
            output(i).formatData=[];
            output(i).type='overlay';
            output(i).defaultDisplayMethod=@(x) FATracksDisplay('Linewidth', 1, 'Color', [255/255 153/255 51/255]); 
            i = ii-7; output(i).name='Focal Adhesion Tracks'; 
            output(i).var='trackFA';
            output(i).formatData=[];
            output(i).type='overlay';
            output(i).defaultDisplayMethod=@(x) FATracksDisplay('Linewidth', 1, 'Color', 'b'); 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Adhesion Boundaries
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
            i = ii-9; output(i).name='Focal Adhesion Boundary';
            output(i).var='adhboundary_FA';
            output(i).formatData=[];
            output(i).type='overlay';
            output(i).defaultDisplayMethod=@(x) AdhBoundaryDisplay('Color', 'b');
            
            i = ii-8; output(i).name='Focal Contact Boundary';
            output(i).var='adhboundary_FC';
            output(i).formatData=[];
            output(i).type='overlay';
            output(i).defaultDisplayMethod=@(x) AdhBoundaryDisplay('Color', [255/255 153/255 51/255]);

        end       

        function funParams = getDefaultParams(owner, varargin)

            % MD.getPackage(MD.getPackageIndex('FocalAdhesionPackage')).outputDirectory_
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            ip.addOptional('iChan', 1:numel(owner.channels_),...
               @(x) all(owner.checkChanNum(x)));
            ip.parse(owner, varargin{:});
            
            % Set default parameters
            funParams.OutputDirectory = [ip.Results.outputDir filesep 'AdhesionAnalysis'];
            funParams.ChannelIndex = ip.Results.iChan;

            funParams.ApplyCellSegMask = true;
            funParams.SegCellMaskProc = []; % Specify Process with cell mask output
            funParams.detectedNAProc = []; % Specify FA detection Process index
            funParams.trackFAProc = []; % Specify FA tracking Process index
            funParams.FAsegProc = []; % Specify FA segmentation Process index           

            funParams.onlyEdge = false; 
            funParams.matchWithFA = true; 
            funParams.minLifetime = 5;  % For tracks
            funParams.reTrack = true;
            funParams.getEdgeRelatedFeatures = true;
            funParams.bandwidthNA = 7;
            funParams.minFALengthMicron = 2;
            
            %% TODO - likely will remove this.
            funParams.backupOldResults = true;           
        end

    end
end
