function LG_showDisplacements_callback
% callback to show Displacements in labelgui2

% for 'update'-button: run this callback with a switch (we need to set the
% buttonDownFcn!). First: check whether the menuItem is still checked (just
% in case), and get figure handle. Then goto showDisplacements with
% figureHandle, where we just replace the axes (subplot(x,y,z,'replace')
% Update-button is placed on the figure with this function.

% get handles
[naviHandles, movieWindowHandles] = LG_getNaviHandles;

% decide to launch or close figure
isChecked = get(naviHandles.LG_navi_menuShowDisplacements,'checked');
if strcmp(isChecked,'on')
    % close figure (will uncheck)
    figureHandle = movieWindowHandles.otherWindows.LG_displacementFigure;
    LG_figureCloseReq(figureHandle);
    
    % uncheck again, just to make sure (in case there's no figure, for
        % example)
        set(naviHandles.LG_navi_menuShowDisplacements,'checked','off')
    
    return
else
    
    % close figure (will uncheck)
    figureHandle = movieWindowHandles.otherWindows.LG_displacementFigure;
    if ~isempty(figureHandle)
    LG_figureCloseReq(figureHandle);
    [naviHandles, movieWindowHandles] = LG_getNaviHandles;
    end
    
    % set checkmark
    set(naviHandles.LG_navi_menuShowDisplacements,'checked','on')
end


% get data
idlist = movieWindowHandles.idlist;
dataProperties = movieWindowHandles.dataProperties;
colorMap = movieWindowHandles.colorMap;

% check navigator for figure position
figurePosition = naviHandles.positions.LG_displacementFigure;

% plot testRatios
[figureHandle,objectHandles] = LG_showDisplacements(idlist,dataProperties,figurePosition,colorMap);

% allow navigation by click
set(objectHandles,'ButtonDownFcn','LG_gotoFrameWithClick');

% set closerequest for figure, remember name of menuItem
set(figureHandle,'Tag','LG_displacementFigure');
set(figureHandle,'UserData','LG_navi_menuShowDisplacements')
set(figureHandle,'CloseRequestFcn','LG_figureCloseReq(gcf)');

% store figureHandle
movieWindowHandles.otherWindows.LG_displacementFigure = figureHandle;

% save movieWindowHandles
guidata(movieWindowHandles.LG_movieWindow,movieWindowHandles);