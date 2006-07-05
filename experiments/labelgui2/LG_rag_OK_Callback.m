function LG_rag_OK_Callback(reAssignHandles,doRecalc)
%LG_rag_OK_Callback is the callback for the two ok-buttons in reAssignGUI
%   The function collects all the information needed by LG_reAssignIdlist

% get handles
[naviHandles,movieWindowHandles] = LG_getNaviHandles;

% get idlist etc
idlist = movieWindowHandles.idlist;
goodTimes = movieWindowHandles.idlistData.goodTimes;
dataProperties = movieWindowHandles.dataProperties;

% current time is not current labelgui time, but time that is displayed on
% RAG
currentTime = reAssignHandles.currentTime;


% read handle-values from GUI
pdValues = LG_rag_getPdValues(reAssignHandles.pdHandles);
goodSpotsIdx = reAssignHandles.goodSpotsIdx;
goodTagIdx = reAssignHandles.goodTagIdx;

% find if assignFuture
assignFuture = get(reAssignHandles.LG_rag_futureFrames_rb,'Value');

% only continue if necessary
if all(pdValues == reAssignHandles.originalPdValues)
    return
end

% Change pdValues to tagIndices. Make sure that size of pdValues remains
% constant
gti = [0;goodTagIdx];
pdValues(:) = [gti(pdValues)];

% check if all tags have been assigned. If not (and no recalc), set recalc
% to 1
if ~all(any(pdValues,2)) && doRecalc == 0
    doRecalc = 1;
end

% if tag2 has been assigned but not tag1, swap
swapRows = pdValues(:,1) == 0 & pdValues(:,2) ~= 0;
pdValues(swapRows,[1,2]) = pdValues(swapRows,[2,1]);

% add spotIdx as first column to pdValues
pdValues = [goodSpotsIdx,pdValues];

% reAssign idlist. Always re-estimate
[idlist,success,dataProperties] =...
    LG_reAssignIdlist(idlist,currentTime, goodTimes, ...
    pdValues, assignFuture, 2-doRecalc,dataProperties);
% if all went well, store idlist
if success == 1
    if ~isempty(dataProperties)
        movieWindowHandles.dataProperties = dataProperties;
        movieWindowHandles.dataPropertiesHasChanged = 1;
        guidata(movieWindowHandles.LG_movieWindow,movieWindowHandles);
    end
    LG_loadIdlist(idlist,0);
end
