classdef AdhesionClassificationProcess < DataProcessingProcess
    
    methods (Access = public)
    
        function obj = AdhesionClassificationProcess(owner, varargin)
    
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
                super_args{2} = AdhesionClassificationProcess.getName;
                super_args{3} = @classifyMovieNascentAdhesions;
                
                if isempty(funParams)
                    funParams = AdhesionClassificationProcess.getDefaultParams(owner,outputDir);
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
            outputList=cell(1,10);
            for i=1:9
                outputList{i}=['class' num2str(i) 'Tracks'];
            end
            outputList{10} = 'adhboundary_Classified';
            outputList{11} = 'selectedGroups';
            outputList{12} = 'iClassesAll'; 

            ip =inputParser;
            ip.addRequired('obj');
            ip.addRequired('iChan', @(x) obj.checkChanNum(x));
            ip.addOptional('iFrame', obj.owner_.nFrames_ ,@(x) obj.checkFrameNum(x));
            ip.addParameter('useCache',true,@islogical);
            ip.addParameter('output', outputList{3}, @(x) all(ismember(x,outputList)));
            ip.parse(obj,iChan,varargin{:})
            output = ip.Results.output;
            varargout = cell(numel(output), 1);
            iFrame = ip.Results.iFrame;
            if ischar(output),output={output}; end
            
            % Data loading
%             s = cached.load(obj.outFilePaths_{5,iChan}, '-useCache', ip.Results.useCache, 'tableTracksNA');
            % Persistent works only for double variable or array. Working
            % around ...
            persistent xCoord yCoord refineFAID stateAll startingFrameExtra endingFrameExtra lastFinishTime lastFinishTimeClass
            if isempty(lastFinishTime)
                lastFinishTime = clock; % assigning current time.. This will be definitely different from obj.finishTime_
            end
            
            if ismember(output,outputList(11))
                iOut=1;
                s = load(obj.outFilePaths_{iOut,iChan}, 'idGroup1Selected', 'idGroup2Selected',...
                    'idGroup3Selected', 'idGroup4Selected', 'idGroup5Selected', 'idGroup6Selected',...
                    'idGroup7Selected', 'idGroup8Selected', 'idGroup9Selected');
                varargout{1}=s;
            elseif ismember(output,outputList(12))
                iOut=4;
                iClasses = cached.load(obj.outFilePaths_{iOut,iChan}, '-useCache', ip.Results.useCache,...
                    'idGroup1','idGroup2','idGroup3','idGroup4','idGroup5','idGroup6','idGroup7','idGroup8','idGroup9');
                
                varargout{1}=iClasses;
            elseif ismember(output,outputList(1:10))
                iAdhProc = obj.owner_.getProcessIndex('AdhesionAnalysisProcess');
                adhAnalProc = obj.owner_.getProcess(iAdhProc);
                if isempty(xCoord) || isempty(yCoord) || isempty(refineFAID) ...
                        || isempty(stateAll) || isempty(startingFrameExtra) || ...
                        isempty(endingFrameExtra) || ~all(adhAnalProc.finishTime_==lastFinishTime)
                    s = load(adhAnalProc.outFilePaths_{1,iChan},'metaTrackData');
                    metaTrackData = s.metaTrackData;
                    fString = ['%0' num2str(floor(log10(metaTrackData.numTracks))+1) '.f'];
                    numStr = @(trackNum) num2str(trackNum,fString);
                    trackFolderPath = metaTrackData.trackFolderPath;
                    if ~exist(trackFolderPath,'dir')
                        adhProcFolder=fileparts(adhAnalProc.outFilePaths_{1,iChan});
                        trackFolderPath = [adhProcFolder filesep 'trackIndividual'];
                        metaTrackData.trackFolderPath = trackFolderPath;
                        disp('Updaing metaTrackData with relocated trackIndividual location...')
                        save(adhAnalProc.outFilePaths_{1,iChan},'metaTrackData');
                    end
                    
                    trackIndPath = @(trackNum) [trackFolderPath filesep 'track' numStr(trackNum) '.mat'];
                    for ii=metaTrackData.numTracks:-1:1
                        curTrackObj = load(trackIndPath(ii),'curTrack');
                        curTrack = curTrackObj.curTrack;
                        % final check on state
                        if iscell(curTrack.state)
                            curTrack.state = strcmp(curTrack.state,'BA')+2*strcmp(curTrack.state,'NA')+...
                                3*strcmp(curTrack.state,'FC')+4*strcmp(curTrack.state,'FA')+...
                                5*strcmp(curTrack.state,'ANA')+6*strcmp(curTrack.state,'Out_of_Band');
                        end
                        tracksNA(ii,1) = curTrack;
                        progressText((metaTrackData.numTracks-ii)/metaTrackData.numTracks,'Loading tracksNA') % Update text
                    end
                    s = struct2table(tracksNA);

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
                    lastFinishTime = adhAnalProc.finishTime_;
                end
                
                %I didn't fully investigate this yet, but refineFAID is
                %sometimes cell array. We need to convert this to matrix
                if iscell(refineFAID)
                    maxFrame = max(cellfun(@length,refineFAID_cell));
                    insuffRows = cellfun(@(x) length(x)<maxFrame,refineFAID_cell);
                    for k=find(insuffRows')
                        refineFAID_cell{k} = [refineFAID_cell{k} ...
                                    NaN(1,maxFrame-length(refineFAID_cell{k}))];
                    end
                    refineFAID = cell2mat(refineFAID_cell);
                end                
                
                fileInfo = dir(obj.outFilePaths_{4,iChan});
                if isempty(lastFinishTimeClass)
                    iClasses = cached.load(obj.outFilePaths_{4,iChan}, '-useCache', ip.Results.useCache,...
                        'idGroup1','idGroup2','idGroup3','idGroup4','idGroup5','idGroup6','idGroup7','idGroup8','idGroup9');
                    lastFinishTimeClass = datevec(fileInfo.datenum);
                elseif all(datevec(fileInfo.datenum)==lastFinishTimeClass)
                    iClasses = cached.load(obj.outFilePaths_{4,iChan}, '-useCache', ip.Results.useCache,...
                        'idGroup1','idGroup2','idGroup3','idGroup4','idGroup5','idGroup6','idGroup7','idGroup8','idGroup9');
                    lastFinishTimeClass = datevec(fileInfo.datenum);
                else
                    iClasses = cached.load(obj.outFilePaths_{4,iChan}, '-useCache', false,...
                        'idGroup1','idGroup2','idGroup3','idGroup4','idGroup5','idGroup6','idGroup7','idGroup8','idGroup9');
                end
                idGroupLabel= 1*iClasses.idGroup1 + ...
                                2*iClasses.idGroup2 + ...
                                3*iClasses.idGroup3 + ...
                                4*iClasses.idGroup4 + ...
                                5*iClasses.idGroup5 + ...
                                6*iClasses.idGroup6 + ...
                                7*iClasses.idGroup7 + ...
                                8*iClasses.idGroup8 + 9*iClasses.idGroup9;
                % Note, could do a stack.

                %% Check struct vs table loading           
    %             if isstruct(s)
    %                 s = s.tableTracksNA;
    %             else
    %                 disp('loaded as table');
    %             end

                nTracks = length(xCoord(:,iFrame));
                number = (1:length(xCoord(:,iFrame)))';
                state = stateAll(:,iFrame);
                iiformat = ['%.' '3' 'd'];

                for iout = 1:numel(output)
                    switch output{iout}
                        case 'class1Tracks'  
                            validState = iClasses.idGroup1;
                        case 'class2Tracks'  
                            validState = iClasses.idGroup2;
                        case 'class3Tracks'  
                            validState = iClasses.idGroup3;
                        case 'class4Tracks'  
                            validState = iClasses.idGroup4;
                        case 'class5Tracks'  
                            validState = iClasses.idGroup5;
                        case 'class6Tracks'  
                            validState = iClasses.idGroup6;
                        case 'class7Tracks'  
                            validState = iClasses.idGroup7;
                        case 'class8Tracks'  
                            validState = iClasses.idGroup8;
                        case 'class9Tracks'  
                            validState = iClasses.idGroup9;
                        case 'adhboundary_Classified'
    %                         validState = (s.refineFAID(:,iFrame));
    %                         validState = cellfun(@(x) ~isempty(x(1,iFrame)) & ~isnan(x(1,iFrame)),validState);
                            validState = state == 4 | state == 3; %strcmp(state, 'FA') | strcmp(state,'FC');
                        otherwise
                            error('Incorrect Output Var type');
                    end   
                    if ~isempty(strfind(output{iout},'Tracks'))

                        vars = {'xCoord', 'yCoord', 'number'};
                        validTracks = validState & startingFrameExtra <= iFrame & endingFrameExtra >= iFrame;                    
                        st = table(xCoord(:,1:iFrame), yCoord(:,1:iFrame), number, ...
                                   'VariableNames', {'xCoord', 'yCoord', 'number'});                    

                        varargout{iout}(nTracks, 1) = struct('xCoord', [], 'yCoord', [], 'number', []);
                        varargout{iout}(validTracks, :) = table2struct(st(validTracks, vars));

                    elseif ~isempty(strfind(output{iout},'adhboundary'))                    

    %                     adhBoundary = cellfun(@(x) x{iFrame}, s{validState, 'adhBoundary'}, 'UniformOutput', false);                         
    %                     varargout{iout} = table2struct(table(adhBoundary, number(validState), ...
    %                                                          'VariableNames',{'adhBoundary', 'number'}));                                 
    %                     MD = obj.owner_;
    %                     iAdhAnalProc=MD.getProcessIndex('AdhesionAnalysisProcess');
    %                     adhAnalProc=MD.getProcess(iAdhAnalProc);
                        p=adhAnalProc.funParams_;
                        labelTifPath = [p.OutputDirectory filesep 'labelTifs'];
                        labelAdhesion = imread(strcat(labelTifPath,'/label',num2str(iFrame,iiformat),'.tif'));
    %                     labelAdhesion = bwlabel(labelAdhesion>0,4);
    %                     maxLabel=max(labelAdhesion(:));
    %                     adhBound = cell(1,maxLabel);
    %                     for ii=1:maxLabel
    %                         curAdhBound = bwboundaries(labelAdhesion==ii,4,'noholes');
    %                         adhBound{ii} = curAdhBound{1}; % strongly assumes each has only one boundary
    %                     end
                        adhBound = bwboundaries(labelAdhesion>0,4,'noholes'); % strongly assumes each has only one boundary

                        validAdhState = refineFAID(validState,iFrame); %cellfun(@(x) x(iFrame),refineFAID(validState));
                        varargout{iout}{1} = adhBound(validAdhState); 
                        varargout{iout}{2} = idGroupLabel(validState); %corresponding classes
                    else
                        varargout{iout} = [];
                    end
                end
            end
        end      
    end

    methods (Static)
        function name = getName()
            name = 'Adhesion Classification';
        end

        function h = GUI()
            h = @adhesionClassificationProcessGUI;
        end
        
        function output = getDrawableOutput()
            numGroups=9;
            colors = distinguishable_colors(numGroups,'k');
            % switching colors between group 6 and 9
            tempColor = colors(6,:);
            colors(6,:) = colors(9,:);
            colors(9,:) = tempColor;
            
            i = 1; output(i).name=['Class ' num2str(i) ' Tracks']; 
            output(i).var=['class' num2str(i) 'Tracks'];
            output(i).formatData=[];
            output(i).type='overlay';
            output(i).defaultDisplayMethod=@(x) FATracksDisplay('Linewidth', 0.8, 'Color', colors(i,:)); 

            i = 2; output(i).name=['Class ' num2str(i) ' Tracks']; 
            output(i).var=['class' num2str(i) 'Tracks'];
            output(i).formatData=[];
            output(i).type='overlay';
            output(i).defaultDisplayMethod=@(x) FATracksDisplay('Linewidth', 0.8, 'Color', colors(i,:)); 
            
            i = 3; output(i).name=['Class ' num2str(i) ' Tracks']; 
            output(i).var=['class' num2str(i) 'Tracks'];
            output(i).formatData=[];
            output(i).type='overlay';
            output(i).defaultDisplayMethod=@(x) FATracksDisplay('Linewidth', 0.8, 'Color', colors(i,:)); 
            
            i = 4; output(i).name=['Class ' num2str(i) ' Tracks']; 
            output(i).var=['class' num2str(i) 'Tracks'];
            output(i).formatData=[];
            output(i).type='overlay';
            output(i).defaultDisplayMethod=@(x) FATracksDisplay('Linewidth', 0.8, 'Color', colors(i,:)); 
            
            i = 5; output(i).name=['Class ' num2str(i) ' Tracks']; 
            output(i).var=['class' num2str(i) 'Tracks'];
            output(i).formatData=[];
            output(i).type='overlay';
            output(i).defaultDisplayMethod=@(x) FATracksDisplay('Linewidth', 0.8, 'Color', colors(i,:)); 
            
            i = 6; output(i).name=['Class ' num2str(i) ' Tracks']; 
            output(i).var=['class' num2str(i) 'Tracks'];
            output(i).formatData=[];
            output(i).type='overlay';
            output(i).defaultDisplayMethod=@(x) FATracksDisplay('Linewidth', 0.8, 'Color', colors(i,:)); 
            
            i = 7; output(i).name=['Class ' num2str(i) ' Tracks']; 
            output(i).var=['class' num2str(i) 'Tracks'];
            output(i).formatData=[];
            output(i).type='overlay';
            output(i).defaultDisplayMethod=@(x) FATracksDisplay('Linewidth', 0.8, 'Color', colors(i,:)); 
            
            i = 8; output(i).name=['Class ' num2str(i) ' Tracks']; 
            output(i).var=['class' num2str(i) 'Tracks'];
            output(i).formatData=[];
            output(i).type='overlay';
            output(i).defaultDisplayMethod=@(x) FATracksDisplay('Linewidth', 0.8, 'Color', colors(i,:)); 
            
            i = 9; output(i).name=['Class ' num2str(i) ' Tracks']; 
            output(i).var=['class' num2str(i) 'Tracks'];
            output(i).formatData=[];
            output(i).type='overlay';
            output(i).defaultDisplayMethod=@(x) FATracksDisplay('Linewidth', 0.8, 'Color', colors(i,:)); 
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Adhesion Boundaries
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
            i = 10; output(i).name='Focal Adhesion Boundary';
            output(i).var='adhboundary_Classified';
            output(i).formatData=[];
            output(i).type='overlay';
            output(i).defaultDisplayMethod=@(x) AdhBoundaryDisplay('Color', 'b');
        end       

        function funParams = getDefaultParams(owner, varargin)

            % MD.getPackage(MD.getPackageIndex('FocalAdhesionPackage')).outputDirectory_
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            adhAnalProc = owner.getProcess(owner.getProcessIndex('AdhesionAnalysisProcess'));
            pAnal=adhAnalProc.funParams_;
            
            ip.addOptional('iChan', pAnal.ChannelIndex,...
               @(x) all(owner.checkChanNum(x)));
            ip.parse(owner, varargin{:});
            
            % Set default parameters
            funParams.OutputDirectory = [ip.Results.outputDir filesep 'AdhesionClassification'];
            funParams.ChannelIndex = ip.Results.iChan;

            funParams.ApplyCellSegMask = true;
            funParams.SegCellMaskProc = []; % Specify Process with cell mask output
            funParams.detectedNAProc = []; % Specify FA detection Process index
            funParams.trackFAProc = []; % Specify FA tracking Process index`````
            funParams.FAsegProc = []; % Specify FA segmentation Process index           

            funParams.labelData=[];
            funParams.useAutomaticallySelectedData = true;
            funParams.manualLabeling=false;
            funParams.useSimpleClassification=false;
            funParams.startingDist=2; % in micron
            %% TODO - likely will remove this.
            funParams.backupOldResults = true;           
        end

    end
end
