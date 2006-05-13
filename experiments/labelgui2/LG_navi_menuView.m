function LG_navi_menuView
% labelgui2-callback executed when the menu File is opened.
% Turns menu items on/off (revert, save)

% get handles
[naviHandles,movieWindowHandles] = LG_getNaviHandles;

% check if data is here
if isempty(movieWindowHandles) || isempty(movieWindowHandles.idlist)
    
    set(naviHandles.LG_navi_menuShowTestRatios,'Enable','off')
    set(naviHandles.LG_navi_menuShowIntensities,'Enable','off')
    set(naviHandles.LG_navi_menuShowDisplacements,'Enable','off')
    set(naviHandles.LG_navi_menuShowDistances,'Enable','off')
    set(naviHandles.LG_navi_menuShowMovieData,'Enable','off')
else
    
    set(naviHandles.LG_navi_menuShowTestRatios,'Enable','on')
    set(naviHandles.LG_navi_menuShowIntensities,'Enable','on')
    set(naviHandles.LG_navi_menuShowDistances,'Enable','on')
    set(naviHandles.LG_navi_menuShowDisplacements,'Enable','on')
    set(naviHandles.LG_navi_menuShowMovieData,'Enable','on')
end
