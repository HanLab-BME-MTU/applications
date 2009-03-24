function fsmDataViewer

% Check fsmDataViewer is not already open
h = findall(0, 'Tag', 'fsmDataViewer');

if ~isempty(h) && ishandle(h)
    set(h, 'Visible', 'on');
    return;
end

% First time 'fsmDataViewer' so call the setup window.
settings = fsmDataViewerSetup;

if (isempty(settings))
    % The user has cancel
    return;
end

% Set up the main window
h = figure('Tag', 'fsmDataViewer',...
    'Name', 'fsmDataViewer',...
    'NumberTitle', 'off');

% Add a slider
sliderStep = [1 5] / settings.numFrames;

uicontrol(h, 'Style', 'slider', ...
    'Units', 'normalized',...
    'Value', sliderStep(1), ...
    'Min', sliderStep(1), ...
    'Max', 1, ...
    'SliderStep', sliderStep, ...
    'Callback', 'sliderShowFrame_Callback', ...
    'Tag', 'sliderShowFrame', ...
    'Position', [0,0,1,0.05]);

% Add menu (set settings, export menu, etc.)

% Attach the settings to the figure
set(h, 'UserData', settings);

% Display the first frame
displayFrame(h, 1);

end