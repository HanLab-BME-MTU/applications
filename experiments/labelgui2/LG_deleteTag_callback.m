function LG_deleteTag_callback
%LG_deleteTag_callback is the callback to delete tags

% get handles
[naviHandles, movieWindowHandles] = LG_getNaviHandles;

% get idlist, goodTimes
idlist = movieWindowHandles.idlist;
goodTimes = movieWindowHandles.idlistData.goodTimes;
currentTime = LG_getCurrentTime;

% find current tag:
% read tag index in userData of current object
tagIdx = get(gco,'UserData');

% ask user again
ans = myQuestdlg(...
    sprintf('Do you really want to remove tag ''%s''?',...
    idlist(1).stats.labelcolor{tagIdx}),'Careful!','Yes','No','Yes');
if ~strmatch(ans,'Yes')
    % if user doesn't agree: return
    return
end

% look at flags. If we're removing a tag that's a good spot, recalc
% completely
flagList = catStruct(1,sprintf('idlist.linklist(%i,5)',tagIdx));
if ~all(ismember(flagList,[2,3]))
    doRecalc = 1;
else
    doRecalc = 0;
end

% delete tag
idlist = LG_deleteTag(idlist,tagIdx,goodTimes);

% recalc if necessary
if doRecalc
    dataProperties = movieWindowHandles.dataProperties;
    % in case the number of spot changes: Ask user about recalc
    if dataProperties.MAXSPOTS ~= nnz(idlist(currentTime).linklist(:,2)>0)
        recalcOpt = [];
    else
        recalcOpt = {0};
    end
    [idlist] = ...
        LG_recalcIdlist(idlist,dataProperties,{0});
end

% save idlist, idlistData - simply reload without replacing everything
LG_loadIdlist(idlist, 0);
