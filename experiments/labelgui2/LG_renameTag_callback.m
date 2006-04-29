function LG_renameTag_callback(newTagName)
%LG_renameTag is the callback to change a tag label

[naviHandles, movieWindowHandles] = LG_getNaviHandles;

% get idlist
idlist = movieWindowHandles.idlist;

% read tag index in userData of current object
tagIdx = get(gco,'UserData');

% change tag name
idlist = LG_renameTag(idlist,tagIdx,newTagName);

% save idlist
LG_loadIdlist(idlist, 0);