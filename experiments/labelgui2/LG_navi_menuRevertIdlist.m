function LG_navi_menuRevertIdlist(hObject,eventdata,naviHandles)
%LG_navi_menuRevertIdlist is the menu callback to revert to the saved idlist

% get handles (yes, I know that I could get them from the input)
[naviHandles, movieWindowHandles] = LG_getNaviHandles;

% load saved idlist
idlist = movieWindowHandles.safeIdlist;

% dataHasChanged:
% 0 : same as original
% 1 : different from original, different from last saved
% 2 : same as original, different from saved
% 3 : different from original, saved
% In principle, it should be impossible to hit revert in cases 0 and 2, but
% who knows how many future errors there will be
% -> if we loaded from the outside, we want dataHasChanged to be at least 1
if naviHandles.launchedFromOutside
    dataHasChanged = 1;
else
    if ~isEven(movieWindowHandles.dataHasChanged)
        dataHasChanged = ...
            movieWindowHandles.dataHasChanged - 1;
    end
end

% revert idlist. LoadIdlist plots
success = LG_loadIdlist(idlist, 0, [], []);

% save dataHasChanged. Reload movieWindowHandles, though. Otherwise, we
% save the old idlist
[naviHandles, movieWindowHandles] = LG_getNaviHandles;
movieWindowHandles.dataHasChanged = dataHasChanged;
guidata(movieWindowHandles.LG_movieWindow, movieWindowHandles);
