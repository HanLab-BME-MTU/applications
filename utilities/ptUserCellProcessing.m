function [newCoord] = ptUserCellProcessing (inputImage, coord)
% ptUserCellProcessing displays an image with coordinates overlayed as red dots. It allows the user to manually edit missed or faulty coordinates.
%
% SYNOPSIS       [newCoord] = ptUserCellProcessing (inputImage, coord)
%
% INPUT          inputImage  : image to process
%                coord       : list of coordinates
%
% OUTPUT         newCoord : edited list of coordinates
%
% DEPENDENCIES   ptUserCellProcessing uses { nothing }
%                                  
%                ptUserCellProcessing is used by { ptInitializeJob }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Colin Glass           Feb 04          Initial release
% Andre Kerstens        Mar 04          Cleaned up source

% Show the image in a figure
initialFigure = figure, imshow (inputImage, []), title ('Left-klick on the unmarked nuclei, right-klick on the dots that are not nuclei // Press ENTER twice when finished.');

% Show the coordinates overlayed as red dots in the same figure
hold on;
plot(coord(:,1),coord(:,2),'r.');

% The user can click on cells. Everytime he or she presses ENTER, the changes will become visible.
% The user can continue to mark cells as long as he or she marks at least one cell between 
% pressing ENTER twice. Then the user will be asked if this is really correct. Left-click indicates yes, 
% and the function returns the new coordinates. Rightclick will let the user continue clicking around.

% Initialize the while loop
found = [];
allAddCoord = [];

% This is a loop that can potentially take forever (just keep on clicking)
while 1
   % Make sure everything takes place on the original figure
   hold on;

   % Get x- and y-coordinates from the mouse (coordinates the user clicked on)
   [xMouseCoord, yMouseCoord, mouseButton] = ginput;

   % Stop if no mouse clicks were provided (and enter was pressed)
   if isempty (xMouseCoord)
       
      % Close initial figure (don't worry we'll draw a new one right away)
      close (initialFigure);
          
      % Remove the cells marked as faulty (right-clicked)
      if ~isempty (found)
         coord (found,:) = [];
      end
                 
      % Add the user-selected extra cells
      %allAddCoord(1,:) = [];
      coord = cat (1, coord, allAddCoord);
   
      % Ask the user if this is really correct and show the final image
      finalFigure = figure, imshow (inputImage, []), title ('Finished? Yes: left-click -- No: right-click, followed by pressing ENTER.');

      % Overlay the red dots again
      hold on;
      plot (coord(:,1), coord(:,2), 'r.');

      % And get input from the user
      [dummy, dummy2, mouseButton] = ginput;

      % Close the current figure window
      hold off;
      
      % mouseButton contains 1 if the user made a left-click. In this case the functions returns
      if ~isempty (mouseButton) & mouseButton(1) == 1
         newCoord = coord;
         close (finalFigure);
         return;
      else
         % We still want to close the figure and continue
         close (finalFigure);
      end
   
      % In case the user didn't left-click, the function continues offering the user to change coordinates. 
      % Since allAddCoord and the coordinates to be removed up to this point are already integrated,
      % we have to reinitialize these two variables
      allAddCoord = [0,0];
      found = [];
   
      % Plot the updated figure on the screen again and ask the user to change coordinates
      initialFigure = figure, imshow (inputImage, []), title ('Left-klick on the unmarked nuclei, right-klick on the dots that are not nuclei // Press ENTER twice when finished.');

      % Also overlay the coordinates in red dots again
      hold on;
      plot (coord (:,1), coord (:,2), 'r.');

      % And get user input
      clear mouseButton;
      clear xMouseCoord;
      clear yMouseCoord;
      [xMouseCoord, yMouseCoord, mouseButton] = ginput;

      hold off;
      %close;
   end
   
   % Continue as long as the user is still clicking around on the image
   left=[];
   right=[];

   % Distinguish left from rightclicks and store the respective coordinates
   left = find (mouseButton == 1);      % 1 means left-clicks
   right = find (mouseButton == 3);     % 3 means right-clicks
       
   % In case left-clicks were made, save the coordinates to be added and plot them
   if ~isempty (left)
      % Get the x and y coordinates to be added
      xAdd = round (xMouseCoord (left));
      yAdd = round (yMouseCoord (left));

      % Make sure we're still plotting on the same figure
      hold on;
      
      % Plot these using a red dot (r.)
      plot (xAdd, yAdd, 'r.');

      % Put the x and y coordinates together in one vector
      coordAdd = cat (2, xAdd, yAdd);

      % Add the newly found coordinate pair to the total
      allAddCoord = cat (1, allAddCoord, coordAdd);
   end

   % In case right-clicks were made, save the coordinates to be removed and plot them
   if ~isempty (right)
      % Get the x and y coordinates to be removed
      xRemove = round (xMouseCoord (right));
      yRemove = round (yMouseCoord (right));
   
      % Search for the actual existing red dots in a confined area
      for index = 1 : length (xRemove)
         [dummy, found(end+1)] = min (abs (coord (:,1) - xRemove (index)) + abs (coord (:,2) - yRemove (index)));
      end
                  
      %coordRemove = [0,0];
      coordRemove = coord (found,:);
      
      % Make sure we're still plotting on the same figure
      hold on;

      % Plot the coordinates to be removed with a yellow dot (y.)
      plot (coordRemove (:,1), coordRemove (:,2), 'y.');

      % Next time a new figure will be drawn
      hold off;
   end

   % Next time a new figure will be drawn
   hold off;
end
