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

end