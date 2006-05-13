function LG_showMovieData_callback(update)
% callback to show testRatios in labelgui2
%
% update: whether to update (1) or to open new figure ({0})


if nargin == 0 || isempty(update)
    update = 0;
end

% get handles
[naviHandles, movieWindowHandles] = LG_getNaviHandles;


% decide to launch or close figure
isChecked = get(naviHandles.LG_navi_menuShowMovieData,'checked');
% check if updating (yes, we could solve this via switch.)
if update
    if strcmp(isChecked,'off')
        % close figure
        figureHandle = movieWindowHandles.otherWindows.LG_movieDataFigure;
        LG_figureCloseReq(figureHandle);

        % uncheck again, just to make sure (in case there's no figure, for
        % example)
        set(naviHandles.LG_navi_menuShowMovieData,'checked','off')

        return
    end
else
    if strcmp(isChecked,'on')
        % close figure (will uncheck)
        figureHandle = movieWindowHandles.otherWindows.LG_movieDataFigure;
        LG_figureCloseReq(figureHandle);

        % uncheck again, just to make sure (in case there's no figure, for
        % example)
        set(naviHandles.LG_navi_menuShowMovieData,'checked','off')

        return
    else

        % close figure (will uncheck)
        figureHandle = movieWindowHandles.otherWindows.LG_movieDataFigure;
        if ~isempty(figureHandle)
            LG_figureCloseReq(figureHandle);
            [naviHandles, movieWindowHandles] = LG_getNaviHandles;
        end

        % set checkmark
        set(naviHandles.LG_navi_menuShowMovieData,'checked','on')
    end
end

%% LOAD DATA

% need: idlist, dataProperties, currentTime, loadedFrames, maxSpots, tagPopupMenuH
idlist = movieWindowHandles.idlist;
dataProperties = movieWindowHandles.dataProperties;
currentTime = LG_getCurrentTime;
if isempty(movieWindowHandles.loadMovieStruct)
    loadedFrames = [1,dataProperties.movieSize(4)];
else
loadedFrames = movieWindowHandles.loadMovieStruct.loadedFrames;
end
maxSpots = max(movieWindowHandles.idlistData.nSpots);




% check navigator for figure position or movieWindow for figure handle
tagPopupMenuH = [];
if update
    figureHandleOrPos = movieWindowHandles.otherWindows.LG_movieDataFigure;
    if ~ishandle(figureHandleOrPos)
        figureHandleOrPos = naviHandles.positions.LG_movieDataFigure;
        update = 0;
    else
        fHandles = guidata(figureHandleOrPos);
        tagPopupMenuH = fHandles.tagPopupMenuH;
    end
else
    figureHandleOrPos = naviHandles.positions.LG_movieDataFigure;
end

% plot testRatios
figureHandle = LG_showMovieData(...
    idlist, dataProperties, currentTime, loadedFrames, maxSpots,...
    figureHandleOrPos, tagPopupMenuH);

if ~update
    % set closerequest for figure, remember name of menuItem
    set(figureHandle,'Tag','LG_testRatiosFigure');
    set(figureHandle,'UserData','LG_navi_menuShowMovieData')
    set(figureHandle,'CloseRequestFcn','LG_figureCloseReq(gcf)');

    % store figureHandle
    movieWindowHandles.otherWindows.LG_movieDataFigure = figureHandle;

    % save movieWindowHandles
    guidata(movieWindowHandles.LG_movieWindow,movieWindowHandles);
end