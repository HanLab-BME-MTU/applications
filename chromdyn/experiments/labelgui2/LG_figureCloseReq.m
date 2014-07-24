function LG_figureCloseReqFcn(hObject,eventdata,handles)
%LG_figureCloseReqFcn is the general closeRequestFunction for labelgui2
% Here, the position gets stored, and the figure deleted

% in case the closeReq is called erroneously: return
if isempty(hObject) || ~ishandle(hObject)
    return
end

% find which figure we're currently closing
currentFigure = get(hObject,'Tag');

% store position, remove handle from movieWindowHandles
[naviHandles,movieWindowHandles] = LG_getNaviHandles(0);
if ~isempty(naviHandles)
    currentFigureHandle = movieWindowHandles.otherWindows.(currentFigure);
    naviHandles.positions.(currentFigure) = get(hObject,'Position');
    guidata(naviHandles.LG_navigator, naviHandles);
    movieWindowHandles.otherWindows.(currentFigure) = [];
    guidata(movieWindowHandles.LG_movieWindow,movieWindowHandles);
    
    % check whether to remove "checked" from menu
    menuHandle = get(hObject,'UserData');
    if ~isempty(menuHandle)
        set(naviHandles.(menuHandle),'checked','off')
    end
else
    % of course, without any other window, this function can only be called
    % by the right gui.
    currentFigureHandle = hObject;
end


% Hint: delete(hObject) closes the figure
delete(currentFigureHandle);