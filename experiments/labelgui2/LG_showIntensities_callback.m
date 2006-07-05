function LG_showIntensities_callback(update)
% callback to show Intensities in labelgui2


% get handles
[naviHandles, movieWindowHandles] = LG_getNaviHandles;

if nargin == 0 || isempty(update)
    update = false;
end

% decide to launch or close figure

isChecked = get(naviHandles.LG_navi_menuShowIntensities,'checked');
if update
    if strcmp(isChecked,'off')
        % close figure
        figureHandle = movieWindowHandles.otherWindows.LG_intensityFigure;
        LG_figureCloseReq(figureHandle);

        % uncheck again, just to make sure (in case there's no figure, for
        % example)
        set(naviHandles.LG_navi_menuShowIntensities,'checked','off')

        return
    end
else
if strcmp(isChecked,'on')
    % close figure (will uncheck)
    figureHandle = movieWindowHandles.otherWindows.LG_intensityFigure;
    LG_figureCloseReq(figureHandle);
    
    % uncheck again, just to make sure (in case there's no figure, for
        % example)
        set(naviHandles.LG_navi_menuShowIntensities,'checked','off')
    
    return
else
    
    % close figure (will uncheck)
    figureHandle = movieWindowHandles.otherWindows.LG_intensityFigure;
    if ~isempty(figureHandle)
    LG_figureCloseReq(figureHandle);
    [naviHandles, movieWindowHandles] = LG_getNaviHandles;
    end
    
    % set checkmark
    set(naviHandles.LG_navi_menuShowIntensities,'checked','on')
end
end


% get data
idlist = movieWindowHandles.idlist;
idlistData = movieWindowHandles.idlistData;
colorMap = movieWindowHandles.colorMap;

% if update, get figure handle, otherwise, get figure position
if update
    figureHandleOrPos = movieWindowHandles.otherWindows.LG_intensityFigure;
    if ~ishandle(figureHandleOrPos)
        figureHandleOrPos = naviHandles.positions.LG_intensityFigure;
        update = 0;
    end
else
    % check navigator for figure position
    figureHandleOrPos = naviHandles.positions.LG_intensityFigure;
end

% plot testRatios
[figureHandle,objectHandles] = LG_showIntensities(idlist,idlistData,colorMap,figureHandleOrPos);

% allow navigation by click
set(objectHandles,'ButtonDownFcn','LG_gotoFrameWithClick');

if ~update
% set closerequest for figure, remember name of menuItem
set(figureHandle,'Tag','LG_intensityFigure');
set(figureHandle,'UserData','LG_navi_menuShowIntensities')
set(figureHandle,'CloseRequestFcn','LG_figureCloseReq(gcf)');
end

% store figureHandle
movieWindowHandles.otherWindows.LG_intensityFigure = figureHandle;

% save movieWindowHandles
guidata(movieWindowHandles.LG_movieWindow,movieWindowHandles);
