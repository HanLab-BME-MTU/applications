function LG_deleteSingleOccurences_callback
%LG_deleteSingleOccurences is the callback for deleting single occurences

% get handles
[naviHandles, movieWindowHandles] = LG_getNaviHandles;

% get idlist, goodTimes
idlist = movieWindowHandles.idlist;
goodTimes = movieWindowHandles.idlistData.goodTimes;


% ask user again
ans = myQuestdlg(...
    'Do you really want to remove all single occurences?',...
    'Careful!','Yes','No','Yes');
if isempty(strmatch(ans,'Yes'))
    % if user doesn't agree: return
    return
end

% delete single occurences. It's a bit silly, but I try to adhere to the
% standards, here
idlist = LG_deleteSingleOccurences(idlist,goodTimes);

% write back idlist
% save idlist, idlistData - simply reload without replacing everything
LG_loadIdlist(idlist, 0);

