function LG_showTestRatios_callback
% callback to show testRatios in labelgui2

% get handles
[naviHandles, movieWindowHandles] = LG_getNaviHandles;

% decide to launch or close figure
isChecked = get(naviHandles.LG_navi_menuShowTestRatios,'checked');
if strcmp(isChecked,'on')
    % close figure (will uncheck)
    figureHandle = movieWindowHandles.otherWindows.LG_testRatiosFigure;
    LG_figureCloseReq(figureHandle);
    
        % uncheck again, just to make sure (in case there's no figure, for
        % example)
        set(naviHandles.LG_navi_menuShowTestRatios,'checked','off')
    
    return
else
    
    % close figure (will uncheck)
    figureHandle = movieWindowHandles.otherWindows.LG_testRatiosFigure;
    if ~isempty(figureHandle)
    LG_figureCloseReq(figureHandle);
    [naviHandles, movieWindowHandles] = LG_getNaviHandles;
    end
    
    % set checkmark
    set(naviHandles.LG_navi_menuShowTestRatios,'checked','on')
end

% get dataProperties
dataProperties = movieWindowHandles.dataProperties;

% get testRatios. Keep it all simple and try to load it from file
% look in idlistDir first
% label_showTestRatios is actually prettier, but there is no time for this.
idlistDir = movieWindowHandles.idlistDir;
try
    load(fullfile(idlistDir,sprintf('testRatios_%s.mat',dataProperties.name)))
catch
    try
        load(fullfile(idlistDir,'testRatios.mat'))
    catch
        % if not working, try movieDir
        movieDir = movieWindowHandles.movieDir;
        try
            load(fullfile(movieDir,'testRatios.mat'))
        catch
            % if this still didn't work, there is no testRatios.
            h = errordlg(sprintf('Sorry, no testRatios found for %s',dataProperties.name));
            uiwait(h);
            % uncheck!
            set(naviHandles.LG_navi_menuShowTestRatios,'checked','off')
            return
        end
    end
end

% check navigator for figure position
figurePosition = naviHandles.positions.LG_testRatiosFigure;

% plot testRatios
[figureHandle,objectHandles] = LG_showTestRatios(testRatios,dataProperties,figurePosition);

% allow navigation by click
set(objectHandles,'ButtonDownFcn','LG_gotoFrameWithClick');

% set closerequest for figure, remember name of menuItem
set(figureHandle,'Tag','LG_testRatiosFigure');
set(figureHandle,'UserData','LG_navi_menuShowTestRatios')
set(figureHandle,'CloseRequestFcn','LG_figureCloseReq(gcf)');

% store figureHandle
movieWindowHandles.otherWindows.LG_testRatiosFigure = figureHandle;

% save movieWindowHandles
guidata(movieWindowHandles.LG_movieWindow,movieWindowHandles);