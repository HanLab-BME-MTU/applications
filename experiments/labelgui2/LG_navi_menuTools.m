function LG_navi_menuTools
% labelgui2-callback executed when the menu File is opened.
% Turns menu items on/off (revert, save)

% get handles
[naviHandles,movieWindowHandles] = LG_getNaviHandles;

% check if data is here
if isempty(movieWindowHandles) || isempty(movieWindowHandles.idlist)
    % disable revert, save, loadIdlist
    set(naviHandles.LG_navi_relinkTags,'Enable','off')
    set(naviHandles.LG_navi_menuDeleteFrame,'Enable','off')
    set(naviHandles.LG_navi_menuDeleteFrameList,'Enable','off')
    set(naviHandles.LG_navi_menuRecalcIdlist,'Enable','off')
    
else
    % disable revert, save, enable load
    set(naviHandles.LG_navi_relinkTags,'Enable','on')
    set(naviHandles.LG_navi_menuDeleteFrame,'Enable','on')
    set(naviHandles.LG_navi_menuDeleteFrameList,'Enable','on')
    set(naviHandles.LG_navi_menuRecalcIdlist,'Enable','on')
    
end
