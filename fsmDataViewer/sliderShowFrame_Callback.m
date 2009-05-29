function sliderShowFrame_Callback

hFig = findall(0, '-regexp','Name','fsmDataViewer');
hSlider = findobj(hFig, 'Tag', 'sliderShowFrame');
settings = get(hFig, 'UserData');

sliderValue = get(hSlider, 'Value');
iFrame = round(sliderValue * settings.numFrames);

% Unlock imtool axes children.
set(hFig, 'HandleVisibility', 'on');

% Update the title of the window.
set(hFig, 'Name', ['fsmDataViewer: frame (' num2str(iFrame) '/'...
    num2str(settings.numFrames) ')' ]);


% Get the axes of hFig.
hAxes = get(hFig, 'CurrentAxes');

% Get the children of axes.
hContent = get(hAxes, 'Children');

% Find the image object in the list of children.
hImage = findobj(hContent, 'type', 'image');

if numel(hImage) ~= 1
    error('Current figure contain none or invalid image data.');
end

% Load and display the new frame
cdata = loadChannels(settings, iFrame);
clim = [min(cdata(:)) max(cdata(:))];
map = get(hFig, 'Colormap');

set(hImage, 'CData', cdata);
set(hFig, 'Colormap', map);
set(hAxes, 'CLim', clim);

% Display layers
displayLayers(hFig, iFrame);

refresh(hFig);

% Relock imtool axes children.
set(hFig, 'HandleVisibility', 'callback');

end