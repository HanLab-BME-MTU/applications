function fsmDataViewer(varargin)

% Check fsmDataViewer is not already open
h = findall(0, '-regexp', 'Name', 'fsmDataViewer');

if ~isempty(h) && ishandle(h)
    set(h, 'Visible', 'on');
    return;
end

% First time 'fsmDataViewer' so call the setup window.
settings = fsmDataViewerSetup(varargin);

if (isempty(settings))
    % The user has cancel
    return;
end

% Set up the main figure (using imtool)
hFig = imtool(loadChannels(settings, 1), []);

% Change 'save as...' menu to custom function
hSaveMenu = findobj(hFig, '-regexp', 'Tag', 'save as menu item');
set(hSaveMenu, 'Callback', @saveMovie);

% Unlock imtool axes children.
set(hFig, 'HandleVisibility', 'on');

% Set the title of the window
set(hFig, 'Name', ['fsmDataViewer: frame (' num2str(1) '/' num2str(settings.numFrames) ')' ]);

% Add a slider
if settings.numFrames > 1
    sliderStep = [1 5] / (settings.numFrames - 1);

    uicontrol(hFig, 'Style', 'slider', ...
        'BackgroundColor', [.91 .91 .91],...
        'Value', 1, ...
        'Min', 1, ...
        'Max', settings.numFrames, ...
        'SliderStep', sliderStep, ...
        'Callback', 'sliderShowFrame_Callback', ...
        'Tag', 'sliderShowFrame', ...
        'Position', [1 40 200, 14]);
end

% Attach the settings to the figure.
set(hFig, 'UserData', settings);
 
% Display layers of the first frame
displayLayers(hFig, 1);

% Relock imtool axes children.
set(hFig, 'HandleVisibility', 'callback');

function saveMovie(varargin)

% Get figure handler
h = findall(0, '-regexp', 'Name', 'fsmDataViewer');
settings = get(h, 'UserData');
nFrames = settings.numFrames;
fString = strcat('%0',num2str(ceil(log10(nFrames)+1)),'.f');

filterSpec = {...
    '*.eps', 'EPS file (*.eps)'; ...
    '*.tif','TIFF image (*.tif)'; ...
    '*.mov','QuickTime movie (*.mov)'};

[fileName,pathName,filterIndex] = uiputfile(filterSpec, 'Save As', 'Untitled.eps');

if ischar(fileName) && ischar(pathName)
    isQuickTime = strcmp(filterSpec{filterIndex,1}, '*.mov');

    if ~isQuickTime
        [path,body] = getFilenameBody(fullfile(pathName,fileName));
        fString = strcat('%0',num2str(ceil(log10(nFrames)+1)),'.f');
        ext = filterSpec{filterIndex,1}(2:end);
        switch filterIndex
            case 1
                driver = '-depsc2';
            case 2
                driver = '-dtiffn';            
        end
    end
    
    [iFirst,iLast]=fsmTrackSelectFramesGUI(1,nFrames,0,'Select frames to be processed:');
    
    scrsz = get(0,'ScreenSize');
    hFig = figure('Position',scrsz);
    set(hFig,'UserData',settings);
    
    if isQuickTime
        evalString = ['MakeQTMovie start ''' fullfile(pathName,fileName) ''''];
        eval(evalString);
    end
    
    if settings.numChannels == 1
        colormap 'gray';
    end
        
    for iFrame = iFirst:iLast  
        % Display channels
        I = loadChannels(settings, iFrame);
        
        imagesc(I,'Parent', gca),axis image,axis off;
        % Display layers
        displayLayers(hFig, iFrame);
        
        if isQuickTime
            MakeQTMovie addaxes
        else
            print(hFig, driver, [path filesep body num2str(iFrame,fString) ext]);
        end
    end
    
    if isQuickTime
         MakeQTMovie finish
    end
    
    close(hFig);
end

