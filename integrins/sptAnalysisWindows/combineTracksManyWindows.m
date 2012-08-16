function [eventCombinedTracks,eventWindowsList,eventNumContributions] = ...
    combineTracksManyWindows(eventCombinedWindows,windowTrackAssignExt)


%Khuloud Jaqaman, July 2011

%get number of activity types
numType = length(eventCombinedWindows);

%get fields in the structure eventCombinedWindows
eventField = fieldnames(eventCombinedWindows(1).eventInfo(1));
numFields = length(eventField);

%get size of window track assignment matrix
[numBands,numSlices,~,numFramesMinus1] = size(windowTrackAssignExt);

%reserve memory for output
for iField = 1 : numFields
    eventCombinedTracks.(eventField{iField}) = {[]};
end
eventCombinedTracks = repmat(eventCombinedTracks,numType,1);
[eventWindowsList,eventNumContributions] = deal(eventCombinedTracks);

%go over activity types
for iType = 1 : numType
    
    %get event information
    eventInfo = eventCombinedWindows(iType).eventInfo;
    
    %go over all fields in eventCombinedWindows (indicating various
    %series relative to activity onset)
    if ~isempty(eventInfo)
        for iField = 1 : numFields
            
            %get maximum increment, whether positive or negative
            windowTrackFrameIndx = eventInfo(1).(eventField{iField});
            [maxInc,numCol] = size(windowTrackFrameIndx);
            
            %reserve memory
            tracksCurrent = cell(maxInc,numCol);
            windowsListCurrent = tracksCurrent;
            numEventsContr = NaN(maxInc,numCol);
            
            %get all relevant windows at all increments
            windowTrackFrameIndxAll = vertcat(eventInfo.(eventField{iField}));
            
            %go over all increments
            for iCol = 1 : numCol
                for iInc = 1 : maxInc
                    
                    %get windows relevant to this increment
                    windowTrackFrameIndx = windowTrackFrameIndxAll(iInc:maxInc:end,iCol);
                    
                    %get number of activity onsets that contribute to this
                    %increment
                    winContributes = zeros(length(windowTrackFrameIndx),1);
                    for iWin = 1 : length(windowTrackFrameIndx)
                        winContributes(iWin) = ~isempty(windowTrackFrameIndx{iWin});
                    end
                    numEventsContr(iInc,iCol) = sum(winContributes);
                    
                    %put window information into an array form
                    windowTrackFrameIndx = vertcat(windowTrackFrameIndx{:});
                    
                    %get tracks in those windows
                    if isempty(windowTrackFrameIndx)
                        tracksCurrent{iInc,iCol} = [];
                        windowsListCurrent{iInc,iCol} = [];
                    else
                        linearInd = sub2ind([numBands numSlices numFramesMinus1 numFramesMinus1],...
                            windowTrackFrameIndx(:,1),windowTrackFrameIndx(:,2),...
                            windowTrackFrameIndx(:,3),windowTrackFrameIndx(:,4));
                        assignmentTmp = windowTrackAssignExt(linearInd);
                        tracksCurrent{iInc,iCol} = vertcat(assignmentTmp{:});
                        windowsListCurrent{iInc,iCol} = windowTrackFrameIndx;
                    end
                    
                end
            end
            
            %store output
            eventCombinedTracks(iType).(eventField{iField}) = tracksCurrent;
            eventWindowsList(iType).(eventField{iField}) = windowsListCurrent;
            eventNumContributions(iType).(eventField{iField}) = numEventsContr;
            
        end %(for iField = 1 : numFields)
    end %(if ~isempty(eventInfo))
    
    %this is a quick hack to fix the number of events count
    if ~iscell(eventNumContributions(iType).onset)
        eventNumContributions(iType).onset(:) = eventNumContributions(iType).onset(1);
        eventNumContributions(iType).befStatic = repmat(eventNumContributions(iType).befStatic(:,1),...
            1,size(eventNumContributions(iType).befStatic,2));
        eventNumContributions(iType).befDynamic = eventNumContributions(iType).befStatic;
        eventNumContributions(iType).aftStatic = repmat(eventNumContributions(iType).aftStatic(:,1),...
            1,size(eventNumContributions(iType).aftStatic,2));
        maxInc = length(eventNumContributions(iType).aftDynamicComb);
        tmp = repmat(eventNumContributions(iType).aftStatic(:,1),1,size(eventNumContributions(iType).aftDynamic,2));
        for iInc = 1 : maxInc
            tmp(iInc,iInc+1:end) = 0;
        end
        eventNumContributions(iType).aftDynamic = tmp;
        eventNumContributions(iType).aftDynamicComb = (sum(tmp))';
    end
    
end %(for iType = 1 : numType)
