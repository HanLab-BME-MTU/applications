function ptShowSlidingFrames
% ptShowSlidingFrames finds the rigth image and coordinates. These it shows in the
% figure opened by ptManualPostProcessJob
%
% SYNOPSIS       ptShowSlidingFrames
%
% INPUT          none (it gets values from the slider created in ptManualPostProcessJob)
%
% OUTPUT         none (it updates the handles object directly)
%
% DEPENDENCIES   ptShowSlidingFrames uses {nothing}
%                                  
%                ptShowSlidingFrames is used by { ptManualPostProcessJob }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Colin Glass           Feb 04          Initial release
% Andre Kerstens        Mar 04          Cleaned up source
% Andre Kerstens        May 04          Renamed to ptShowSlidingFrames.m


% this is the callback of the slider, created in manualpostpro
% What we do here is:
% - find out which frame the user currently wants to look at
% - show this frame in the figure (created in manualpostpro)
% - plot the coordinates of the cells into this picture
% - plot the cell numbers (near the respective cells) with the 
%   callback manrelink
% In case the user clicks on the number of the cell, manrelink comes into play.

% Delete the axes currently shown in the figure on the screen
delete (gca)

% Look for objects needed for information
sliderHandle = findall (0, 'Tag', 'pictureslide');
%hObject = findall (0, 'Tag', 'GUI_add_pb');
hObject = findall (0, 'Tag', 'GUI_filelist_lb');
frameCounterHandle = findall (0, 'Style', 'text', 'Tag', 'picturecount');

% Use the hObject just found to get to the handles structure
handles = guidata (hObject);

% Fetch the jobvalues and image directory
imageDirectory = handles.jobData(1).imagefilepath;
imageName      = handles.jobData(1).imagename;
firstImage     = handles.jobData(1).firstimg;
lastImage      = handles.jobData(1).lastimg;
increment      = handles.jobData(1).increment;
imageRange     = handles.ma;
imageNameList  = handles.jobData(1).imagenameslist;
intensityMax   = handles.jobData(1).intensitymax;

% Get the current value of the slider, so that we know which frame the user wants to process
sliderValue = get (sliderHandle, 'Value');
sliderValue = round (sliderValue * imageRange);

% Calculate the frame number to show
imageNumber = (sliderValue - 1) * increment + firstImage;

% Write the current frame number in the little window above the slider
set (frameCounterHandle, 'String', num2str (imageNumber));

% Read the image frame from disk
cd (imageDirectory);
fileName = char (imageNameList (imageNumber));
image = imreadnd2 (fileName, 0, intensityMax);

% Show the frame on the screen in the current figure
hold on;
imshow (image, []), title (num2str (imageNumber));
hold off;

% Get the cells corresponding to this frame. We create a third column
% (first two being [x,y]) with the row indices of MPM. We are NOT interested in
% the indices the cells will have in the vector cellsWithNums!!! Why?
% because in MPM every cell has it's own row, so the row indices of MPM is
% the actual number of the cell

% Identify the real cells (at least one coord different from zero)
realCellIndex = find (handles.allMPM{1}(:, 2 * sliderValue - 1) | handles.allMPM{1}(:, 2 * sliderValue));

% Find the row indices from a transposed MPM matrix
cellsWithNums = zeros (size (handles.allMPM{1}, 1), 3);
cellsWithNums(:,3) = [1:1:size(handles.allMPM{1},1)]';

% Grab all rows in MPM, so that the row indices correspond to the cells
cellsWithNums(:,1:2) = handles.allMPM{1}(:, 2 * sliderValue - 1 : 2 * sliderValue);

% Now take the cells identified as real cells (at least one coord different from zero)
% and plot those as red dots. The cell number is written as colored text on the current axes.
hold on;
plot (cellsWithNums (realCellIndex, 1), cellsWithNums (realCellIndex, 2), 'r.');
txt = text (cellsWithNums (realCellIndex, 1), cellsWithNums (realCellIndex, 2), num2str (cellsWithNums (realCellIndex, 3)), 'Color', 'r');
hold off;

% Depending on who was the caller, set some object values
if handles.whichcallback == 1
   set (txt, 'ButtonDownFcn', 'manrelink'); 
elseif  handles.whichcallback == 2
   if ~isempty (handles.selectedcells)
      handles.selectedcells = [];
      % Update the handles structure
      guidata (hObject, handles);
   end
   set(txt,'ButtonDownFcn','cellselect');
end

% That's it: wait for the next user action
