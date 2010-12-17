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
set(hFig, 'Name', ['fsmDataViewer: frame (' num2str(1) '/' num2str(settings.nFrames) ')' ]);

% Add a slider
if settings.nFrames > 1
    sliderStep = [1 5] / (settings.nFrames - 1);

    uicontrol(hFig, 'Style', 'slider', ...
        'BackgroundColor', [.91 .91 .91],...
        'Value', 1, ...
        'Min', 1, ...
        'Max', settings.nFrames, ...
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


filterSpec = {'*.tif','TIFF image (*.tif)'};

[fileName,pathName] = uiputfile(filterSpec, 'Save As', 'Untitled.tif');

if ischar(fileName) && ischar(pathName)
    % Get figure handler
    h = findall(0, '-regexp', 'Name', 'fsmDataViewer');
    settings = get(h, 'UserData');
    nFrames = settings.nFrames;

    % Get the crop area if any
    aspectRatio = [1 1];
    hCrop = findobj(get(h, 'CurrentAxes'), 'Tag', 'imrect');
    if isempty(hCrop) || ~ishandle(hCrop)
        crop = false;
    else
        crop = true;
        hC = get(hCrop, 'Children');
        
        x = arrayfun(@(h) get(h, 'XData'), hC(1:8));
        y = arrayfun(@(h) get(h, 'YData'), hC(1:8));
        coords = [x, y];
        cropCenter = mean(coords, 1);
        cropSize = max(coords) - min(coords);
        if max(coords) == coords(1)
            aspectRatio = [max(cropSize) / min(cropSize) 1];
        else
            aspectRatio = [1 max(cropSize) / min(cropSize)];
        end
    end
    
    % Create the figure where stack will be painted onto.
    hFig = figure(...
        'Visible', 'off',...
        'Renderer', 'painters',...
        'Position', [100 100 ceil(aspectRatio * 1000)],...
        'PaperPositionMode', 'auto',...
        'UserData', settings);
    hAxes = axes;
    set(hAxes, 'Position', [0 0 1 1]);
    set(hFig,'CurrentAxes',hAxes);
    
    [path,body] = getFilenameBody(fullfile(pathName,fileName));
    fString = strcat('%0',num2str(ceil(log10(nFrames)+1)),'.f');
    tmpOutputFile = fullfile(path, 'tmp.eps');
    
    [iFirst,iLast]=fsmTrackSelectFramesGUI(1,nFrames,0,'Select frames to be processed:');
    
    if numel(settings.channels) == 1
        colormap 'gray';
    end
    
    for iFrame = iFirst:iLast
        cla(hAxes);
        
        % Display channels
        I = loadChannels(settings, iFrame);
        imagesc(I,'Parent', hAxes);
        
        % Display layers
        displayLayers(hFig, iFrame);

        axis(hAxes,'off','image');

        if crop
            % Change the view point of the camera
            set(hAxes, 'CameraPositionMode', 'manual');
            set(hAxes, 'CameraTargetMode', 'manual');
            set(hAxes, 'CameraUpVectorMode', 'manual');
            set(hAxes, 'CameraViewAngleMode', 'manual');

            pos = get(hAxes, 'CameraPosition');
            theta = get(hAxes, 'CameraViewAngle') * pi / 180;
            pos(3) = max(cropSize) / (2 * tan(theta/2));
            pos(1:2) = cropCenter;

            set(hAxes, 'CameraPosition', pos);
            set(hAxes, 'CameraTarget', [pos(1:2), 0]);
        end
        
        % Write an .eps file
        print(hFig, '-depsc2', '-r300', '-noui', tmpOutputFile);
        % Convert the .eps into .tiff using imageMagick
        eval(['!convert -colorspace rgb ' tmpOutputFile ' tiff:' fullfile(path, [body, num2str(iFrame,fString), '.tif'])]);
    end
    
    % Remove the tmp.eps
    if exist(tmpOutputFile, 'file');
        delete(tmpOutputFile);
    end
    
    close(hFig);
end

