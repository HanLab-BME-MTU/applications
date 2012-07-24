function eventCombinedTracks = combineTracksManyWindows(...
    eventCombinedWindows,windowTrackAssignExt)


%Khuloud Jaqaman, July 2011

%get number of activity types
numType = length(eventCombinedTracks);

%go over activity types
for iType = 1 : numType
    
    %get event information
    eventInfo = eventCombinedWindows.eventInfo;
    numEvent = size(eventInfo,1);
    
    %tracks at protrusion onset
    windowsAssignIndx = vertcat(eventInfo.indxWindowAssignOnset);
    tracksAll = 
    
    
end
