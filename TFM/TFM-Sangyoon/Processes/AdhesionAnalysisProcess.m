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
        
        function output = loadChannelOutput(obj, iChan, varargin)
            % Input check
            outputList = {'detBA','detectedFA','detFA','detFC','detNA',...
                          'adhboundary','tracks','staticTracks'};

            ip =inputParser;
            ip.addRequired('obj');
            ip.addRequired('iChan', @(x) obj.checkChanNum(x));
            ip.addOptional('iFrame', [] ,@(x) obj.checkFrameNum(x));
            ip.addParameter('useCache',true,@islogical);
            ip.addParameter('output', outputList{1}, @(x) all(ismember(x,outputList)));
            ip.parse(obj,iChan,varargin{:})
            output = ip.Results.output;
            iFrame = ip.Results.iFrame;
            if ischar(output),output={output}; end
            
            % Data loading
            % load outFilePaths_{1,iChan} (should contain tracksNA)
            s = cached.load(obj.outFilePaths_{1,iChan}, '-useCache', ip.Results.useCache, 'tableTracksNA');
            if isstruct(s)
                s = s.tableTracksNA;
            else
                disp('loaded as table');
            end
            %% Check struct vs table loading
            xCoord = s.xCoord(:,iFrame);
            yCoord = s.yCoord(:,iFrame);
            pres = s.presence(:,iFrame);
            state = categorical(s.state(:,iFrame));

            t = table(xCoord, yCoord, state, pres);

            for output_sel = output
                switch output_sel{1}
                    case 'detectedFA'  
                        output = t;
                    case 'detBA' 
                        rows = t.state == 'BA';
                        vars = {'xCoord','yCoord'};
                        output = t{rows,vars};                 
                    case 'detNA'
                        rows = t.state == 'NA';
                        vars = {'xCoord','yCoord'};
                        output = t{rows,vars};                 
                    case 'detFC'
                        rows = t.state == 'FC';
                        vars = {'xCoord','yCoord'};
                        output = t{rows,vars};                 
                    case 'detFA'
                        rows = t.state == 'FA';
                        vars = {'xCoord','yCoord'};
                        output = t{rows,vars};                 
                    case 'adhboundary'
                    case 'tracks'
                    case 'staticTracks'
                    otherwise
                        error('Incorrect Output Var type');
                end   
            end
        end      
    end

% t = tableTracksNA;
% for i = 1:50
%     xCoord = t.xCoord(:,i);
%     yCoord = t.yCoord(:,i);
%     pres = t.presence(:,i);
%     state = categorical(t.state(:,i));
%     iFrame = repmat(i,[length(xCoord),1]);
%     ti = table(xCoord, yCoord, state, pres, iFrame);
%     if i == 1
%         tb = ti;
%     else
%         tb = vertcat(tb,ti);
%     end
  
% end
% -------------> need stack!

    methods (Static)
        function name = getName()
            name = 'Focal Adhesion Analysis';
        end

        function h = GUI()
            h = @focalAdhesionAnalysisProcessGUI;
        end
        
        function output = getDrawableOutput()
            % output(1).name = 'Adhesion Detections';
            % keySet = {'BA','NA','FC','FA'};
            % valueSet = {'g','r', 'y', 'b'};
            % ColorsDict = containers.Map(keySet,valueSet);
            % % Detection Points
            % % Colors = 'grob'; % optional
            % output(1).var = 'detectedFA';
            % output(1).formatData = [];%obj.formatSpotOutput;
            % output(1).type = 'overlay';
            % output(1).defaultDisplayMethod=@(x)FASpotDisplay('ColorDict', ColorsDict);
            
            % BA Detection (Before Adhesion)
            % colors = 'g';
            i = 1;
            output(i).name='Before Adhesion Detection'; 
            output(i).var='detBA';
            output(1).formatData=[];
            output(1).type='overlay';
            output(1).defaultDisplayMethod=@(x) LineDisplay('Marker','d',...
                'LineStyle','none', 'Color', 'g');            
            
            % NA Detection Nascent Adhesion
            i = 2;
            output(i).name='Nascent Adhesion Detection'; 
            output(i).var='detNA';
            output(i).formatData=[];
            output(i).type='overlay';
            output(i).defaultDisplayMethod=@(x) LineDisplay('Marker','o',...
                'LineStyle','none', 'LineWidth', .7, 'Color', 'r'); 
            
            % FC Detection Focal Contact
            i = 3;
            output(i).name='Focal Contact Detection'; 
            output(i).var='detFC';
            output(i).formatData=[];
            output(i).type='overlay';
            output(i).defaultDisplayMethod=@(x) LineDisplay('Marker','o',...
                'LineStyle','none', 'LineWidth', .75, 'Color', [255/255 153/255 51/255]); 

            % FA Detection Focal Adhesion
            i = 4;
            output(i).name='Focal Adhesion Detection'; 
            output(i).var='detFA';
            output(i).formatData=[];
            output(i).type='overlay';
            output(i).defaultDisplayMethod=@(x) LineDisplay('Marker','o',...
                'LineStyle','none', 'LineWidth', .75, 'Color', 'b'); 

            % BA Track
            % NA Track
            % FC Track
            % FA Track


            % % Adhesion Boundaries
            % colorsAdhBound = hsv(numel(obj.owner_.channels_));
            % output(2).name='Adh. Boundaries';
            % output(2).var='adhboundary';
            % output(2).formatData=@DetectionProcess.formatOutput;
            % output(2).type='overlay';
            % output(2).defaultDisplayMethod=@(x) LineDisplay('Marker','o',...
            %     'LineStyle','none','Color', colorsAdhBound(x,:));
            
            % % Tracks
            % colors = hsv(numel(obj.owner_.channels_));
            % output(3).name='FA Tracks';
            % output(3).var='tracks';
            % output(3).formatData=@TrackingProcess.formatTracks;
            % output(3).type='overlay';
            % output(3).defaultDisplayMethod = @(x)FATracksDisplay(...
            %     'Color',colors(x,:));
            
            % output(4).name='Static FA tracks';
            % output(4).var='staticTracks';
            % output(4).formatData=@TrackingProcess.formatTracks;
            % output(4).type='overlay';
            % output(4).defaultDisplayMethod = @(x)FATracksDisplay(...
            %     'Color',colors(x,:), 'useDragtail', false);



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
            
            %% TODO - likely will remove this.
            funParams.backupOldResults = true;           
        end

        function displayTracks = formatTracks(tracks)
            % Format tracks structure into compound tracks for display
            % purposes
            
            % Determine the number of compound tracks
            nCompoundTracks = cellfun('size',{tracks.tracksCoordAmpCG},1)';
            trackIdx = 1:length(tracks);
            
            % Filter by the tracks that are nonzero
            filter = nCompoundTracks > 0;
            tracks = tracks(filter);
            trackIdx = trackIdx(filter);
            nCompoundTracks = nCompoundTracks(filter);
            
            % Get the track lengths (nFrames x 8)
            trackLengths = cellfun('size',{tracks.tracksCoordAmpCG},2)';
            % Unique track lengths for batch processing later
            uTrackLengths = unique(trackLengths);
            % Running total of displayTracks for indexing
            nTracksTot = [0 cumsum(nCompoundTracks(:))'];
            % Total number of tracks
            nTracks = nTracksTot(end);
            
            % Fail fast if no track
            if nTracks == 0
                displayTracks = struct.empty(1,0);
                return
            end
            
            % Number of events in each seqOfEvents for indexing
            % Each will normally have 2 events, beginning and end
            nEvents = cellfun('size',{tracks.seqOfEvents},1);

            % Initialize displayTracks structure
            % xCoord: x-coordinate of simple track
            % yCoord: y-coordinate of simple track
            % events (deprecated): split or merge
            % number: number corresponding to the original input compound
            % track number
            % splitEvents: times when splits occur
            % mergeEvents: times when merges occur
            displayTracks(nTracks,1) = struct('xCoord', [], 'yCoord', [], 'number', [], 'splitEvents', [], 'mergeEvents', []);
            
            hasSeqOfEvents = isfield(tracks,'seqOfEvents');
            hasLabels = isfield(tracks, 'label');
            
            if(hasLabels)
                labels = vertcat(tracks.label);
                if(size(labels,1) == nTracks)
                    hasPerSegmentLabels = true;
                end
            end

            % Batch by unique trackLengths
            for trackLength = uTrackLengths'
                %% Select original track numbers
                selection = trackLengths == trackLength;
                sTracks = tracks(selection);
                sTrackCoordAmpCG = vertcat(sTracks.tracksCoordAmpCG);
                % track number relative to input struct array
                sTrackIdx = trackIdx(selection);
                % index of selected tracks
                siTracks = nTracksTot(selection);
                
                % runs in current selection
                snCompoundTracks = nCompoundTracks(selection);
                snTracksTot = [0 ; cumsum(snCompoundTracks)];
                
                % decode run lengths
                snTracks = snTracksTot(end);
                idx = zeros(snTracks,1);
                idx(snTracksTot(1:end-1) + 1) = ones(size(snCompoundTracks));
                idx = cumsum(idx);

                % absolute index of tracks
                iTracks = ones(snTracks,1);
                iTracks(snTracksTot(1:end-1) + 1) = [siTracks(1); diff(siTracks)'] - [0 ; snCompoundTracks(1:end-1)-1];
                iTracks = cumsum(iTracks) + 1;
                
                % grab x and y coordinate matrices
                xCoords = sTrackCoordAmpCG(:,1:8:end);
                yCoords = sTrackCoordAmpCG(:,2:8:end);
                
                %% Process sequence of events
                if(hasSeqOfEvents)
                    % make sequence of events matrix
                    seqOfEvents = vertcat(sTracks.seqOfEvents);
                    nSelectedEvents = nEvents(selection);
                    iStartEvents = [1 cumsum(nSelectedEvents(1:end-1))+1];

                    % The fifth column is the start frame for each track
                    seqOfEvents(iStartEvents,5) = [seqOfEvents(1) ; diff(seqOfEvents(iStartEvents,1))];
                    seqOfEvents(:,5) = cumsum(seqOfEvents(:,5));

                    % The sixth column is the offset for the current selected
                    % tracks
                    seqOfEvents(iStartEvents,6) = [0; snCompoundTracks(1:end-1)];
                    seqOfEvents(:,6) = cumsum(seqOfEvents(:,6));

                    % Isolate merges and splits
                    seqOfEvents = seqOfEvents(~isnan(seqOfEvents(:,4)),:);

                    % Apply offset 
                    seqOfEvents(:,3) = seqOfEvents(:,3) + seqOfEvents(:,6);
                    seqOfEvents(:,4) = seqOfEvents(:,4) + seqOfEvents(:,6);

                    % Number of Frames
                    nFrames = trackLength/8;

                    %% Splits
                    % The 2nd column indicates split (1) or merge(2)
                    splitEvents = seqOfEvents(seqOfEvents(:,2) == 1,:);
                    % Evaluate time relative to start of track
                    splitEventTimes = splitEvents(:,1) - splitEvents(:,5);
                    % Time should not exceed the number of coordinates we have
                    splitEvents = splitEvents(splitEventTimes < nFrames,:);
                    splitEventTimes = splitEventTimes(splitEventTimes < nFrames,:);
                    iTrack1 = splitEvents(:,3);
                    iTrack2 = splitEvents(:,4);

                    % Use accumarray to gather the splitEventTimes into cell
                    % arrays
                    if(~isempty(splitEventTimes))
                        splitEventTimeCell = accumarray([iTrack1 ; iTrack2], [splitEventTimes ; splitEventTimes ],[size(xCoords,1) 1],@(x) {x'},{});
                    else
                        splitEventTimeCell = cell(size(xCoords,1),1);
                    end

                    leftIdx = (splitEventTimes-1)*size(xCoords,1) + iTrack1;
                    rightIdx = (splitEventTimes-1)*size(xCoords,1) + iTrack2;
                    xCoords(leftIdx) = xCoords(rightIdx);
                    yCoords(leftIdx) = yCoords(rightIdx);

                    %% Merges
                    % The 2nd column indicates split (1) or merge(2)
                    mergeEvents = seqOfEvents(seqOfEvents(:,2) == 2,:);
                    mergeEventTimes = mergeEvents(:,1) - mergeEvents(:,5);
                    mergeEvents = mergeEvents(mergeEventTimes < nFrames,:);
                    mergeEventTimes = mergeEventTimes(mergeEventTimes < nFrames,:);
                    mergeEventTimes = mergeEventTimes + 1;
                    iTrack1 = mergeEvents(:,3);
                    iTrack2 = mergeEvents(:,4);

                    if(~isempty(mergeEventTimes))
                        mergeEventTimeCell = accumarray([iTrack1 ; iTrack2], [mergeEventTimes ; mergeEventTimes ],[size(xCoords,1) 1],@(x) {x'},{});
                    else
                        mergeEventTimeCell = cell(size(xCoords,1),1);
                    end


                    leftIdx = (mergeEventTimes-1)*size(xCoords,1) + iTrack1;
                    rightIdx = (mergeEventTimes-1)*size(xCoords,1) + iTrack2;
                    xCoords(leftIdx) = xCoords(rightIdx);
                    yCoords(leftIdx) = yCoords(rightIdx);
                end
                          
                %% Load cells into struct fields
                for i=1:length(iTracks)
                    iTrack = iTracks(i);
                    displayTracks(iTrack).xCoord = xCoords(i,:);
                    displayTracks(iTrack).yCoord = yCoords(i,:);
                    displayTracks(iTrack).number = sTrackIdx(idx(i));
                end            
                
                if(hasSeqOfEvents)
                    [displayTracks(iTracks).splitEvents] = splitEventTimeCell{:};
                    [displayTracks(iTracks).mergeEvents] = mergeEventTimeCell{:};
                end
                
                if hasLabels
                    if hasPerSegmentLabels
                        for i=1:length(iTracks)
                            iTrack = iTracks(i);
                            displayTracks(iTrack).label = labels(iTrack);
                        end    
                    else
                        [displayTracks(iTracks).label] = sTracks(idx).label;
                    end
                end
            end
        end      

    end
end
