function success = LG_launchMovieWindow(movie, movieDir, dataProperties,loadMovieStruct)
%LG_launchMovieWindow launches a new movie window (but doesn't plot)

success = 0;

% get the necessary handles
[naviHandles, movieWindowHandles] = LG_getNaviHandles;

% find if there's already another window open
if ~isempty(movieWindowHandles)
    % currently, we just close the current movieWindow. Later, we can allow
    % several windows open at the same time
    % (via pulldown in navi-window. Don't forget option 'all')
    windowHasBeenClosed = close(movieWindowHandles.LG_movieWindow);
    if ~windowHasBeenClosed
        % success remains 0. Die silently.
        return
    end
end

% make new figure. We will set the handleVisibility to callback after
% creating the axes
movieWindowH = figure('Name',dataProperties.name,'NumberTitle','off',...
    'Visible','off','Tag','LG_movieWindow',...
    'CloseRequestFcn','LG_movieWindowCloseReq',...
    'MenuBar','figure');

% set position if we haven't stored it yet

if isempty(naviHandles.positions.LG_movieWindow)

    screenSize = get(0,'MonitorPositions');
    % get labelgui- and movieWindow-position in pixels
    naviUnits = get(naviHandles.LG_navigator,'Units');
    movieWindowUnits = get(movieWindowH,'Units');
    set(movieWindowH,'Units','Pixels');
    set(naviHandles.LG_navigator,'Units','Pixels');
    movieWindowPosition = get(movieWindowH,'Position');
    naviPosition = get(naviHandles.LG_navigator,'Position');

    % x = navigator + 7, y = screenHeight - 60, probably because the
    % name-bar is not included in figure-size calculations.
    movieWindowPosition(1) = 7 + naviPosition(1) + naviPosition(3);
    movieWindowPosition(2) = screenSize(4) - 75 - movieWindowPosition(4);
    set(movieWindowH,'Position',movieWindowPosition);

    % reset pixels
    set(movieWindowH,'Units',movieWindowUnits);
    set(naviHandles.LG_navigator,'Units',naviUnits);

else
    set(movieWindowH,'Position',naviHandles.positions.LG_movieWindow)
end


% --- create axes --- 

% to make the rectangle-zoom: take this part into a new function, start by
% clearing the children of the current figure!

% calculate frame size in microns to use for imshow later.
frameSizeMu = [dataProperties.PIXELSIZE_XY, dataProperties.PIXELSIZE_XY,...
    dataProperties.PIXELSIZE_Z] .*...
    dataProperties.movieSize(1:3);
% the two axes only cover 0.85 of the figure
% relAxesSize(1,1:2) = frameSizeMu([1,3])/sum(frameSizeMu([1,3])) * 0.94;
% relAxesSize(2,1:2) = frameSizeMu([2,3])/sum(frameSizeMu([2,3])) * 0.94;
% relAxesSize = frameSizeMu./sum(frameSizeMu) * 0.94
[dummy,mxyi]=max(frameSizeMu(1:2));
relAxesSize = frameSizeMu/sum(frameSizeMu([mxyi,3])) * 0.94;

% Set aspect ratios, and initialize image so that we later just replace the
% image values.
xyAxesH = ...
    axes('Position',[0.02, 0.04+relAxesSize(3),...
    relAxesSize(2), relAxesSize(1)]);
set(xyAxesH,'PlotBoxAspectRatio',[relAxesSize(2),relAxesSize(1),1],...
    'NextPlot','add','Visible','off');
xyImageH = imagesc([0,frameSizeMu(2)]+0.5*dataProperties.PIXELSIZE_XY,...
    [0,frameSizeMu(1)]+0.5*dataProperties.PIXELSIZE_XY,...
    rand(dataProperties.movieSize([1,2])));
axis image

yzAxesH = ...
    axes('Position',[0.04+relAxesSize(2), 0.04+relAxesSize(3), ...
    relAxesSize(3), relAxesSize(1)],...
    'NextPlot','add','Visible','off');
yzImageH = imagesc([0,frameSizeMu(3)]+0.5*dataProperties.PIXELSIZE_Z,...
    [0,frameSizeMu(1)]+0.5*dataProperties.PIXELSIZE_XY,...
    rand(dataProperties.movieSize([1,3])));
axis image

set(yzAxesH,'PlotBoxAspectRatio',[relAxesSize(3),relAxesSize(1),1]);
xzAxesH = ...
    axes('Position',[0.02, 0.02, relAxesSize(2),relAxesSize(3)]);
set(xzAxesH,'PlotBoxAspectRatio',[relAxesSize(2),relAxesSize(3),1],...
    'NextPlot','add','Visible','off');
xzImageH = imagesc([0,frameSizeMu(2)]+0.5*dataProperties.PIXELSIZE_XY,...
    [0,frameSizeMu(3)]+0.5*dataProperties.PIXELSIZE_Z,...
    rand(dataProperties.movieSize([3,2])));
axis image

colormap('gray')



% now that the axes have been created, we make the figure handle invisible
set(movieWindowH,'HandleVisibility','callback');

% store all the info

% movieWindow. Store data and initialize all fields
movieWindowHandles = guidata(movieWindowH);
movieWindowHandles.xyAxesH = xyAxesH;
movieWindowHandles.yzAxesH = yzAxesH;
movieWindowHandles.xzAxesH = xzAxesH;
movieWindowHandles.xyImageH = xyImageH;
movieWindowHandles.yzImageH = yzImageH;
movieWindowHandles.xzImageH = xzImageH;
movieWindowHandles.frameSizeMu = frameSizeMu;
movieWindowHandles.movie = movie;
movieWindowHandles.movieDir = movieDir;
movieWindowHandles.loadMovieStruct = loadMovieStruct;
movieWindowHandles.dataProperties = dataProperties;
movieWindowHandles.dataHasChanged = []; % initialize only
movieWindowHandles.dataPropertiesHasChanged = 0;
movieWindowHandles.idlist = [];
movieWindowHandles.idname = [];
movieWindowHandles.idlistDir = [];
movieWindowHandles.dataFileName = [];
movieWindowHandles.colorMap = [];
movieWindowHandles.idlistData = [];
movieWindowHandles.flagData = [];

% store handles of other windows. To make life easier in the
% closeRequestFcn, store movieWindowH twice. 
% !! if you're adding a field here, add it in labelgui2-OpeningFcn, too
% Important: Assign movieWindow last!!
movieWindowHandles.otherWindows.LG_reAssignGUI = [];
movieWindowHandles.otherWindows.LG_movieDataFigure = [];
movieWindowHandles.otherWindows.LG_testRatiosFigure = [];
movieWindowHandles.otherWindows.LG_intensityFigure = [];
movieWindowHandles.otherWindows.LG_distanceFigure = [];
movieWindowHandles.otherWindows.LG_displacementFigure = [];
movieWindowHandles.otherWindows.LG_movieWindow = movieWindowH;
% store window handle - a figure does not usually know itself
movieWindowHandles.LG_movieWindow = movieWindowH;

guidata(movieWindowH, movieWindowHandles);

% store movieWindowHandle in navigator
naviHandles.movieWindowH = movieWindowH;
guidata(naviHandles.LG_navigator,naviHandles);

% return success
success = 1;