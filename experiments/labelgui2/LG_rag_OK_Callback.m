function LG_rag_OK_Callback(reAssignHandles,doRecalc)
%LG_rag_OK_Callback is the callback for the two ok-buttons in reAssignGUI
%   The function collects all the information needed by LG_reAssignIdlist

% get handles
[naviHandles,movieWindowHandles] = LG_getNaviHandles;

% get idlist etc
idlist = movieWindowHandles.idlist;
currentTime = LG_getCurrentTime;
goodTimes = movieWindowHandles.idlistData.goodTimes;
dataProperties = movieWindowHandles.dataProperties;

% read handle-values from GUI
pdValues = LG_rag_getPdValues(reAssignHandles.pdHandles);
goodSpotsIdx = reAssignHandles.goodSpotsIdx;

% find if assignFuture
assignFuture = get(reAssignHandles.LG_rag_futureFrames_rb,'Value');

% only continue if necessary
if all(pdValues == reAssignHandles.originalPdValues)
    return
end

% check if all tags are assigned. Return if not. I haven't found a good way
% to deal with both tag deletions and reassignments at the same time
originalValues = unique(reAssignHandles.originalPdValues);
currentValues = unique(pdValues);
if ~all(ismember(originalValues,currentValues))
    h = errordlg(...
        'This tool cannot be used to remove tags. Please assign them all',...
        'Error');
    uiwait(h);
    return
end

% Change pdValues to tagIndices
gsi = [0;goodSpotsIdx];
pdValues = [gsi(pdValues)];

% if tag2 has been assigned but not tag1, swap
swapRows = pdValues(:,1) == 0 & pdValues(:,2) ~= 0;
pdValues(swapRows,[1,2]) = pdValues(swapRows,[2,1]);

% add spotIdx as first column to pdValues
pdValues = [goodSpotsIdx,pdValues];

% reAssign idlist
[idlist,success,dataProperties] =...
    LG_reAssignIdlist(idlist,currentTime, goodTimes, ...
    pdValues, assignFuture, doRecalc,dataProperties);
% if all went well, store idlist
if success == 1
    if ~isempty(dataProperties)
        movieWindowHandles.dataProperties = dataProperties;
        movieWindowHandles.dataPropertiesHasChanged = 1;
        guidata(movieWindowHandles.LG_movieWindow,movieWindowHandles);
    end
    LG_loadIdlist(idlist,0);
end
