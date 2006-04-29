function LG_deleteSpot_callback
%LG_deleteSpot_callback is the callback that assembles all the information necessary to delete a spot

% get handles
[naviHandles, movieWindowHandles] = LG_getNaviHandles;

% get idlist, goodTimes
idlist = movieWindowHandles.idlist;
currentTime = LG_getCurrentTime;
goodTimes = movieWindowHandles.idlistData.goodTimes;
dataProperties = movieWindowHandles.dataProperties;

% find current tag:
% read tag index in userData of current object
tagIdx = get(gco,'UserData');

% find spot to delete
spotNumber = idlist(currentTime).linklist(tagIdx,2);

if spotNumber == 0
    h = errordlg(['This spot can''t be deleted ',...
        'because it doesn''t exist anymore (this is a bug)'],'Bug found');
    uiwait(h);
    return
end

% delete spot. We need dataProperties and goodTimes to call recalc or
% deleteTag.
% Use recalc(3) to estimate position and amplitude if good spot has been
% deleted; use deleteTag to remove single occurence
idlist = LG_deleteSpot(idlist,currentTime,spotNumber,dataProperties,goodTimes);


% save idlist, idlistData - simply reload without replacing everything
LG_loadIdlist(idlist, 0);
