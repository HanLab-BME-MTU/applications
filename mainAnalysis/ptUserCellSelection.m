function [selectedCells] = ptUserCellSelection (frameNr, handles)
% ptUserCellSelection displays an image with coordinates overlayed as red dots. 
% It allows the user to manually select separate cells.
%
% SYNOPSIS       [selectedCells] = ptUserCellSelection (frameNr, handles)
%
% INPUT          frameNr  : image to select from
%                handles  : gui handles struct
%
% OUTPUT         selectedCells : list of cell numbers (not coordinates)
%
% DEPENDENCIES   ptUserCellSelection uses { nothing }
%                                  
%                ptUserCellSelection is used by { Polytrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Feb 05          Initial release

% Load image
cd (handles.jobData(1).imagefilepath);
fileName = char(handles.jobData(1).imagenameslist(frameNr)); 
inputImage = imreadnd2(fileName,0,handles.jobData(1).intensitymax);

% Get already selected cells before
if isfield(handles.jobData(1),'selectedcells')
    if ~isempty(handles.jobData(1).selectedcells)
        prevSelectedCellsIndx = handles.jobData(1).selectedcells(:,frameNr);
        prevSelectedCellsIndx = prevSelectedCellsIndx(find(prevSelectedCellsIndx)); 
    end
end

% Show the image in a figure
selectFig = figure; imshow (inputImage, []);
title ('Select marked nuclei by left-clicking on them with the mouse // Right-click when finished.');

% Get the coordinates from the MPM
cellCoord(:,:) = handles.allMPM{1}(:, (2*frameNr-1):(2*frameNr));
cellCoordAll = cellCoord;
cellCoord = cellCoord(find (cellCoord(:,1) & cellCoord(:,2)),:); 

% Get already selected cells from MPM
if exist('prevSelectedCellsIndx','var')
   prevSelectedCells = handles.allMPM{1}(:, (2*frameNr-1):(2*frameNr));
   %prevSelectedCells(find(prevSelectedCells(:,1) == 0 & prevSelectedCells(:,2) == 0),:) = [];
   prevSelectedCells = prevSelectedCells(prevSelectedCellsIndx,:);
end

% Show the coordinates overlayed as red dots in the same figure
hold on;
plot(cellCoord(:,1),cellCoord(:,2),'r.');

% Plot already selected cells as yellow dots
if exist('prevSelectedCells','var')
   plot(prevSelectedCells(:,1),prevSelectedCells(:,2),'y.');
end

% The user can left-click on nuclei. Everytime the user presses ENTER, the selections will become visible.
% The user can continue to select cells as long as she selects at least one
% cell between clicking and pressing ENTER. A right-click will end the
% process and return the selected cells.

% Initialize the while loop
selectedCell = [];
allSelectedCells = [];

% Loop which will be ended by pressing enter twice
while 1
   % Make sure everything takes place on the original figure
   hold on;

   % Get x- and y-coordinates from the mouse (coordinates the user clicked on)
   [xMouseCoord, yMouseCoord, mouseButton] = ginput;

   % Stop if the right button was clicked (and enter was pressed)
   if ~isempty (mouseButton) & mouseButton(end) == 3
                    
      selectedCells = allSelectedCells;

      % Sort the vector
      selectedCells = sort(selectedCells);

      % Only keep unique entries
      selectedCells = unique(selectedCells);

      % Close the figure
      close (selectFig);
      return;
   end
   
   % Continue as long as the user is still clicking around on the image
   left=[];

   % Store the coordinates of the pixels the user clicked on
   left = find (mouseButton == 1);      % 1 means left-clicks
       
   % Save the coordinates to be added and plot them
   if ~isempty (left)
       
      % Get the x and y coordinates
      xAdd = round(xMouseCoord(left));
      yAdd = round(yMouseCoord(left));

      % Find the nearest real coordinates
      found = [];
      for index = 1 : length(xAdd)
         [dummy,found(end+1)] = min(abs(cellCoordAll(:,1) - xAdd(index)) + abs(cellCoordAll(:,2) - yAdd(index)));
      end
      
      % Transpose the coord matrix
      found = found';
      
      % Make sure we're still plotting on the same figure
      hold on;
      
      % Plot these using a yellow dot (y.)
      if ~isempty(found)
         plot (cellCoordAll(found,1),cellCoordAll(found,2), 'y.');

         % Add the newly found coordinate pair to the total
         allSelectedCells = (cat(1, allSelectedCells, found));
      end
   end

   % Next time a new figure will be drawn
   hold off;
end
