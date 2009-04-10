function sliderShowFrame_Callback

hFig = findall(0, 'Tag', 'fsmDataViewer');
hSlider = findobj(hFig, 'Tag', 'sliderShowFrame');
settings = get(hFig, 'UserData');

sliderValue = get(hSlider, 'Value');
iFrame = round(sliderValue * settings.numFrames);

% Get background of the first frame so that we can create the
% current axes.
B = getBackground(settings, iFrame);

% Unlock imtool axes children.
set(hFig, 'HandleVisibility', 'on');

% Update the title of the window.
set(hFig, 'Name', ['fsmDataViewer: frame (' num2str(iFrame) '/' num2str(settings.numFrames) ')' ]);

% Get the axes of imtool.
hAxes = get(hFig, 'CurrentAxes');

% Get the children of axes.
hContent = get(hAxes, 'Children');

% Find the image object in the list of children.
hImage = findobj(hContent, 'type', 'image');

if numel(hImage) ~= 1
    error('Current figure contain none or invalid image data.');
end

% Update the image data
set(hImage, 'CData', B);

% Relock imtool axes children.
set(hFig, 'HandleVisibility', 'callback');

end