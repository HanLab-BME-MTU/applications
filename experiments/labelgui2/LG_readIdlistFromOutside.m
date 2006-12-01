function [idlist,dataProperties] = LG_readIdlistFromOutside
%LG_READIDLISTFROMOUTSIDE reads the currently saved idlist from labelgui2

% get handles
[naviHandles, movieWindowHandles] = LG_getNaviHandles(0);

if isempty(naviHandles)
    idlist = [];
    dataProperties = [];
else
    % read idlist from movieWindowHandles and return
    idlist = movieWindowHandles.idlist;
    dataProperties = movieWindowHandles.dataProperties;
end

% that's it already.