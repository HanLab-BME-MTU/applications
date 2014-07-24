function LG_deleteFrame_callback
%LG_deleteFrame is the callback for deleting a frame

% get handles
[naviHandles, movieWindowHandles] = LG_getNaviHandles;

% get idlist
idlist = movieWindowHandles.idlist;
dataProperties = movieWindowHandles.dataProperties;

% find currentTime, goodTimes
currentTime = LG_getCurrentTime;
goodTimes = movieWindowHandles.idlistData.goodTimes;

% if idlist(currentTime) has been deleted already, do nothing
if isempty(idlist(currentTime).linklist)
    return
end

% ask whether the user really wants to delete the frame
ans = myQuestdlg(...
    sprintf('Do you really want to delete frame %i?',currentTime),...
    '','Yes & recalc','Yes','No','Yes');

if isempty(ans) || strcmp(ans,'No')
    return
end
if strcmp(ans,'Yes')
    recalc = 0;
else
    recalc = 1;
end

% delete frame
idlist = LG_deleteFrame(idlist,dataProperties,currentTime,goodTimes,recalc);


% save idlist, idlistData - simply reload without replacing everything
LG_loadIdlist(idlist, 0);

        