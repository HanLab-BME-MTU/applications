function LG_recalcIdlist_callback
%LG_recalcIdlist_callback loads the necessary data and calls LG_recalcIdlist

% get handles
[naviHandles, movieWindowHandles] = LG_getNaviHandles;

% get idlist, dataProperties
idlist = movieWindowHandles.idlist;
dataProperties = movieWindowHandles.dataProperties;

% recalc without options (have user choose)
[idlist,dataProperties,success] = ...
    LG_recalcIdlist(idlist,dataProperties);

% if no success: die
if ~success
    return
end

% if dataProperties are returned, update them
if ~isempty(dataProperties)
    movieWindowHandles.dataProperties = dataProperties;
    movieWindowHandles.dataPropertiesHasChanged = 1;
    guidata(movieWindowHandles.LG_movieWindow,movieWindowHandles);
end

% save idlist, idlistData - simply reload without replacing everything
LG_loadIdlist(idlist, 0);