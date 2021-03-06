function ptManualPostProcessJob (hObject)
% ptManualPostProcessJob opens a figure and puts a slider into it; the callback belonging
% to that slider (ptShowSlidingFrames) takes care of showing the correct frame
%
% SYNOPSIS       ptManualPostProcessJob (hObject)
%
% INPUT          hObject : handle of an object of the GUI calling ptManualPostProcessJob
%
% OUTPUT         none
%
% DEPENDENCIES   ptManualPostProcessJob uses {nothing}
%                                  
%                ptManualPostProcessJob is used by { PolyTrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Colin Glass           Feb 04          Initial release
% Andre Kerstens        Mar 04          Cleaned up source


% A routine that, together with it's subprograms, allows the user to
% manually postprocess the results of the main analysis. 
% This programm only opens a figure and puts a slider into it. 
% The sliders callback is ptShowSlidingFrames, so look there if you want to
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

% Fetch the jobvalues
imageDirectory = handles.jobData(1).imagefilepath;
imageName      = handles.jobData(1).imagename;
firstImage     = handles.jobData(1).firstimg;
realFirstImage = handles.jobData(1).realfirstimg;
lastImage      = handles.jobData(1).lastimg;
realLastImage  = handles.jobData(1).reallastimg;
increment      = handles.jobData(1).increment;
imageNameList  = handles.jobData(1).imagenameslist;
intensityMax   = handles.jobData(1).intensitymax;

% Check that we actually have some images to show. If not, give an error
% message
checkName = [imageDirectory filesep char(imageNameList(realFirstImage))];
if exist (checkName, 'file') ~= 2
    errorStr = ['No images available in directory ' imageDirectory];
    h = errordlg(errorStr);
    uiwait(h);          % Wait until the user presses the OK button  
    return;
end

% Calculate the image range taking into account the increment between frames
%imageRange = floor ((lastImage - firstImage + 1) / increment + 0.001);
%imageRange = floor ((lastImage - firstImage + 1) / increment);
imageRange = size (handles.allValidFrames{1},2);
handles.ma = imageRange;

% Generate the slider step values for the uicontrol
% First the arrow slide step (1):
%slider_step(1) = 1 / (imageRange - 1);
slider_step(1) = 1 / imageRange;
% Then the trough step size (5):
%slider_step(2) = 5 / (imageRange - 1);
slider_step(2) = 5 / imageRange;

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
% the function ptShowSlidingFrames when moved
sliderHandle = uicontrol ('Style', 'slider', ...
                          'Units', 'normalized', ... 
                          'Value', 1/(imageRange), ...
                          'Min', 1/(imageRange), ...
                          'Max', 1, ...
                          'SliderStep', slider_step, ...
                          'Callback', 'ptShowSlidingFrames', ...
                          'Tag', 'pictureslide', ...
                          'Position', [0.02,0.02,0.05,0.9]);

% Now we show the first image including red dots for the cells overlaid on the image
% Read the image frame from disk
cd (imageDirectory);
fileName = char (imageNameList (realFirstImage));
image = imreadnd2 (fileName, 0, intensityMax);

% Show the frame on the screen in the current figure
hold on;
%imshow (image, []), title (num2str (firstImage));
imshow (image, [min(image(:)) max(median(image(:)))*2]), title(num2str(firstImage));
hold off;

% Make sure the axis are correct
axis([1 size(image,2) 1 size(image,1)]);

% Get the cells corresponding to the first frame. We create a third column
% (first two being [x,y]) with the row indices of MPM. We are NOT interested in
% the indices the cells will have in the vector cellsWithNums!!! Why?
% because in MPM every cell has it's own row, so the row indices of MPM is
% the actual number of the cell

% Identify the real cells (at least one coord different from zero)
realCellIndex = find (handles.allMPM{1}(:, 1) | handles.allMPM{1}(:, 2));

% Find the row indices from a transposed MPM matrix
cellsWithNums = zeros (size (handles.allMPM{1}, 1), 3);
cellsWithNums(:,3) = [1:1:size(handles.allMPM{1},1)]';

% Grab all rows in MPM, so that the row indices correspond to the cells
cellsWithNums(:,1:2) = handles.allMPM{1}(:,1:2);

% Now take the cells identified as real cells (at least one coord different from zero)
% and plot those as red dots. The cell number is written as colored text on the current axes.
hold on;
dots = plot (cellsWithNums (realCellIndex, 1), cellsWithNums (realCellIndex, 2), 'r.');
set(dots,'Tag','dots','ButtonDownFcn','ptActionOnCells');
txt = text (cellsWithNums (realCellIndex, 1), cellsWithNums (realCellIndex, 2), num2str (cellsWithNums (realCellIndex, 3)), 'Color', 'r');

% That's it: wait for the next user action
hold off;
