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

% Load channels
settings = loadSequence(settings);

imtool(settings.sequence(:, :, :, 1), []);

% Set up the main figure (using imtool)
% hFig = imtool(S(:, :, 1), []);
% set(hFig, 'Tag', 'fsmDataViewer',...
%     'NumberTitle', 'off');
% 
% % Unlock imtool axes children.
% set(hFig, 'HandleVisibility', 'on');
% 
% % Set the title of the window
% set(hFig, 'Name', ['fsmDataViewer: frame (' num2str(1) '/' num2str(settings.numFrames) ')' ]);
% 
% % Add a slider
% sliderStep = [1 5] / settings.numFrames;
% 
% uicontrol(hFig, 'Style', 'slider', ...
%     'Units', 'normalized',...
%     'Value', sliderStep(1), ...
%     'Min', sliderStep(1), ...
%     'Max', 1, ...
%     'SliderStep', sliderStep, ...
%     'Callback', 'sliderShowFrame_Callback', ...
%     'Tag', 'sliderShowFrame', ...
%     'Position', [0,0,1,0.05]);
% 
% % Attach the settings to the figure.
% set(hFig, 'UserData', settings);
% 
% % Display layers of the first frame
% displayLayers(hFig, 1);
% 
% % Relock imtool axes children.
% set(hFig, 'HandleVisibility', 'callback');

end