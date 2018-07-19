function LG_movieWindowCloseReq
%LG_movieWindowCloseReq is the closerequest function for the LG_movieWindow

[naviHandles, movieWindowHandles] = LG_getNaviHandles;

% if naviHandles don't exist, delete figure.
if isempty(naviHandles)
    delete(gcbo)
    return
end

% if figure officially doesn't exist anymore - delete it
if isempty(movieWindowHandles)
    delete(gcbo)
    return
end

% check whether dataHasChanged. If yes, ask user to save. Otherwise, close
% and remove windowH from navigator-guidata
% dataHasChanged:
% 0 : same as original
% 1 : different from original, different from last saved
% 2 : same as original, different from saved
% 3 : different from original, saved
if any(ismember(movieWindowHandles.dataHasChanged,[1,2]))
    ans = myQuestdlg(...
        'Idlist has changed from last saved: Do you want to save it before closing?',...
        '','Yes','No','Cancel','Yes');

    switch ans
        case 'Yes'
            % LG_saveIdlist
        case 'No'
            % do nothing
        case {'','Cancel'}
            % don't close anything: return
            return
    end
end

% remember window positions and delete. This requires that the
% field for the movieWindow handle is created last in LG_launchMovieWindow
windowNames = fieldnames(movieWindowHandles.otherWindows);
for w = 1: length(windowNames)
    windowHandle = movieWindowHandles.otherWindows.(windowNames{w});
    if ishandle(windowHandle)
        % remember position
        naviHandles.positions.(windowNames{w}) = ...
            get(windowHandle,'Position');
        % delete window
        delete(windowHandle);
    end
end

% remove entry in navigator
naviHandles.movieWindowH = [];
% reset title
set(naviHandles.LG_navi_movieName_txt,'String','no movie loaded');
set(naviHandles.LG_navi_flagName_pd,'String','no idlist loaded','Value',1);

% save naviHandles
guidata(naviHandles.LG_navigator,naviHandles);