function sliderShowFrame_Callback

h = findall(0, 'Tag', 'fsmDataViewer');
hSlider = findobj(h, 'Tag', 'sliderShowFrame');
settings = get(h, 'UserData');

sliderValue = get(hSlider, 'Value');
iFrame = round(sliderValue * settings.numImages);

displayFrame(h, iFrame);

end