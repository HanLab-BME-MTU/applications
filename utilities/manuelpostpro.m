function manuelpostpro (hObject)
% manuelpostpro opens a figure and puts a slider into it; the callback belonging
% to that slider (changeframe) takes care of showing the correct frame
%
% SYNOPSIS       manuelpostpro (hObject)
%
% INPUT          hObject : handle of an object of the GUI calling manuelpostpro
%
% OUTPUT         none
%
% DEPENDENCIES   manuelpostpro uses {nothing}
%                                  
%                manuelpostpro is used by { PolyTrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Colin Glass           Feb 04          Initial release
% Andre Kerstens        Mar 04          Cleaned up source


% A routine that, together with it's subprograms, allows the user to
% manually postprocess the results of the main analysis. 
% This programm only opens a figure and puts a slider into it. 
% The sliders callback is changeframe, so better look there if you want to
% know more.

% if the slider (and less important, the little windos with the number of
% the current frame) already exist, delete them now.
imageCounterHandle = findall (0, 'Style', 'text', 'Tag', 'picturecount');
sliderHandle = findall (0, 'Style', 'slider', 'Tag', 'pictureslide');

if ~isempty (imageCounterHandle)
    delete (imageCounterHandle);
end
if ~isempty (sliderHandle)
    delete (sliderHandle);
end

% Get the handles structure
handles = guidata(hObject);

% Initialize the list of cells that will be linked manually (... means empty list)
handles.listofcells='...';

% Fetch the jobvalues and image directory
imageDirectory = handles.postpro.imagepath;
imageName      = handles.jobvalues.imagename;
firstImage     = handles.jobvalues.firstimage;
lastImage      = handles.jobvalues.lastimage;
increment      = handles.jobvalues.increment;

% Calculate the image range taking into account the increment between frames
imageRange = floor ((lastImage - firstImage) / increment + 0.001);
handles.ma = imageRange;

% Generate the slider step values for the uicontrol
slider_step(1) = 1 / (imageRange - 1);
slider_step(2) = 3 / (imageRange - 1);

% Update the handles structure
guidata(hObject, handles);

% Draw a new figure on the screen
figure,

% Draw the frame counter in the figure; it is identified by the tag picturecount
imageCounterHandle = uicontrol ('Style', 'text',...
                                'Units', 'normalized',...
                                'Tag', 'picturecount',...
                                'Position', [0.02,0.93,0.05,0.06]);

% Set the frame counter to the first image number
set (imageCounterHandle, 'String', num2str(firstImage));

% Draw the slider in the figure; it is identified by the tag pictureslide and calls
% the function changeframe when moved
sliderHandle = uicontrol ('Style', 'slider', ...
                          'Units', 'normalized', ... 
                          'Value', 0.00000001, ...
                          'Min', 0, ...
                          'Max', 1, ...
                          'SliderStep', slider_step, ...
                          'Callback', 'changeframe', ...
                          'Tag', 'pictureslide', ...
                          'Position', [0.02,0.02,0.05,0.9]);
