function LG_navi_menuFile
% labelgui2-callback executed when the menu File is opened.
% Turns menu items on/off (revert, save)

% get handles
[naviHandles,movieWindowHandles] = LG_getNaviHandles;

% check if data is here
if isempty(movieWindowHandles) 
    % disable revert, save, loadIdlist
    set(naviHandles.LG_navi_menuRevertIdlist,'Enable','off')
    set(naviHandles.LG_navi_menuSaveIdlist,'Enable','off')
    set(naviHandles.LG_navi_menuLoadIdlist,'Enable','off')
elseif isempty(movieWindowHandles.idlist)
    % disable revert, save, enable load
    set(naviHandles.LG_navi_menuRevertIdlist,'Enable','off')
    set(naviHandles.LG_navi_menuSaveIdlist,'Enable','off')
    set(naviHandles.LG_navi_menuLoadIdlist,'Enable','on')
else
    % allow loading idlist
    set(naviHandles.LG_navi_menuLoadIdlist,'Enable','on')
    % check if dataHasChanged
    % dataHasChanged:
    % 0 : same as original
    % 1 : different from original, different from last saved
    % 2 : same as original, different from saved
    % 3 : different from original, saved
    switch movieWindowHandles.dataHasChanged
        case 0
            set(naviHandles.LG_navi_menuRevertIdlist,'Enable','off')
            set(naviHandles.LG_navi_menuSaveIdlist,'Enable','off')
        case 1
            set(naviHandles.LG_navi_menuRevertIdlist,'Enable','on')
            set(naviHandles.LG_navi_menuSaveIdlist,'Enable','on')
        case 2
            set(naviHandles.LG_navi_menuRevertIdlist,'Enable','off')
            set(naviHandles.LG_navi_menuSaveIdlist,'Enable','on')
        case 3
            set(naviHandles.LG_navi_menuRevertIdlist,'Enable','on')
            set(naviHandles.LG_navi_menuSaveIdlist,'Enable','off')
        otherwise
            h = errordlg('Unknown option for dataHasChanged has been encountered in LG_navi_menuFile');
            uiwait(h)
    end
end
