function ptCalculateNeighbourhoodInteractions (ptPostpro, MPM)
% ptCalculateNeighbourhoodInteractions plots speed information gathered in MPM. 
%
% SYNOPSIS       ptCalculateNeighbourhoodInteractions (ptPostpro, MPM)
%
% INPUT          ptPostpro : a structure which contains the information
%                            from the GUI
%                MPM       : matrix containing the cell tracks
%                
% OUTPUT         None (plots are directly shown on the screen) 
%
% DEPENDENCIES   ptCalculateNeighbourhoodInteractions  uses {nothing}
%                                  
%                ptCalculateNeighbourhoodInteractions is used by { PolyTrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Jul 04          Initial version

% First assign all the postpro fields to a meaningfull variable
startFrame = ptPostpro.firstimg;
endFrame = ptPostpro.lastimg;
increment = ptPostpro.increment;
plotStartFrame = ptPostpro.plotfirstimg;
plotEndFrame = ptPostpro.plotlastimg;
savePath = ptPostpro.saveallpath;
jobPath = ptPostpro.jobpath;
imageName = ptPostpro.imagename;
increment = ptPostpro.increment;
numberOfFrames = ceil((plotEndFrame - plotStartFrame) / increment) + 1;
cellProps = ptPostpro.cellProps;
clusterProps = ptPostpro.clusterProps;
binSize = ptPostpro.binsize;
frameInterval = round (ptPostpro.timeperframe / 60);    % In minutes
pixelLength = ptPostpro.mmpixel;
maxDistance = ptPostpro.maxdistance;
maxDistToNeighbour = maxDistance;

% Over how many frames to we calculate trajectories
% eg 5 frames means 4 trajectories
nrOfTrajectories = ptPostpro.nrtrajectories - 1;

% Initialize the trajectory matrix
avgTrajFrame = zeros (1, numberOfFrames-nrOfTrajectories-1);

% Initialize MPM counter (start with 0 if plotStartFrame=startFrame)
MPMCount = ceil ((plotStartFrame - startFrame) / increment);

% Initialize X-axis vector
xAxis = zeros (1, numberOfFrames-nrOfTrajectories-1);
iCount = 0;

% These calculations can take a while so set the mouse pointer to busy
set(gcf,'Pointer','watch');

% Go through every frame of the set.
for frameCount = plotStartFrame : increment : plotEndFrame
    
   % Increase MPM counter
   MPMCount = MPMCount + 1;
   
   % There should be enough frames left at the end to calc. traj.
   if (MPMCount + nrOfTrajectories) < numberOfFrames
   
      % Store the frame number for display on the x-axis
      iCount = iCount + 1;
      xAxis (iCount) = frameCount;
   
      % Initialize the average trajectory for the frame
      trajFrame = 0;
      
      % Get the cell list for this frame and throw the zero entries out
      cells = MPM (:, 2*MPMCount-1 : 2*MPMCount);
      cells = cells (find (cells (:,1) ~= 0 & cells (:,2) ~= 0),:);
      
      % Triangulate all the cells with their neighbours
      triangleIndex = delaunay (cells(:,1), cells(:,2));
      
      % Initialize counter for cells with nearest-enough neighbours
      cellsWithNeighbours = 0;
      
      % Find neighbours for all these cells and do some trajectory calculations
      for cCount = 1 : size(cells,1)
          
         % Find the entries in triangleIndex for cell 'iCount'
         cellEntries = triangleIndex(find (triangleIndex(:,1)==cCount | triangleIndex(:,2)==cCount | ...
                                          triangleIndex(:,3)==cCount),:);
         
         % Extract the neighbours (incl duplicates)
         neighbours = [];
         nCount = 1;
         for hCount = 1 : size (cellEntries, 1)
            if cellEntries(hCount,1) == cCount
               neighbours(nCount,1) = cellEntries(hCount,2);
               nCount = nCount + 1;
               neighbours(nCount,1) = cellEntries(hCount,3);
               nCount = nCount + 1;
            end
            if cellEntries(hCount,2) == cCount
               neighbours(nCount,1) = cellEntries(hCount,1);
               nCount = nCount + 1;
               neighbours(nCount,1) = cellEntries(hCount,3);
               nCount = nCount + 1;
            end
            if cellEntries(hCount,3) == cCount
               neighbours(nCount,1) = cellEntries(hCount,1);
               nCount = nCount + 1;
               neighbours(nCount,1) = cellEntries(hCount,2);
               nCount = nCount + 1;
            end
         end
         % Take only the unique neighbours
         neighbours = unique (neighbours);
         
         % Calculate the distance to all of the neighbours
         distance = [];
         for dCount = 1 : length (neighbours)   
            distance(dCount,1) = sqrt ((cells(cCount, 1) - cells(neighbours(dCount,1), 1))^2 + ...
                                       (cells(cCount, 2) - cells(neighbours(dCount,1), 2))^2);  
         end
         
         % Only keep the neighbours that are close enough
         neighbours = neighbours (find (distance < maxDistToNeighbour));                                 
                                          
         % In case the cell has any neighbours do some work
         if ~isempty (neighbours)
             
            cellsWithNeighbours = cellsWithNeighbours + 1;
            
            % For every neighbour we are going to calculate avg trajectory
            % length for the specified amount of frames
            trajNeighbours = 0;
            for nCount = 1 : length (neighbours)
               
               % Fetch the neighbour from cells
               curNeighbour = cells (neighbours(nCount,1),:);
               
               % Find the current neighbour in MPM
               mpmIndex = find (MPM(:,2*MPMCount-1) == curNeighbour(:,1) & ...
                                MPM(:,2*MPMCount) == curNeighbour(:,2));
                            
               % Calculate trajectories for the specified amount of frames
               % for this neighbour
               traj = [];
               for tCount = 0 : nrOfTrajectories
                  traj(tCount+1) = sqrt((MPM(mpmIndex, 2*MPMCount-1+(tCount*2)) - MPM(mpmIndex, 2*MPMCount-1+(tCount*2+2)))^2 + ...
                                      (MPM(mpmIndex, 2*MPMCount+(tCount*2)) - MPM(mpmIndex, 2*MPMCount+(tCount*2+2)))^2);
               end
               
               % Average the trajectories and calculate in micrometers
               avgTraj = (sum (traj) / (nrOfTrajectories+1));
               
               % Sum the average traj. per neighbour
               trajNeighbours = trajNeighbours + avgTraj;
            end
            
            % Average the trajectory for all neighbours
            avgTrajNeighbours = trajNeighbours / length (neighbours); 
             
            % Sum the avg trajectory for all neighbours
            trajFrame = trajFrame + avgTrajNeighbours;
             
         end  % if ~isempty (neighbours)
               
      end  % for cCount = 1 : size(cells,1)
      
      % Average the trajectory for all cells with neighbours
      avgTrajFrame(iCount) = (trajFrame / cellsWithNeighbours) * pixelLength;
      
   end   % if (MPMCount + nrOfTrajectories) < numberOfFrames
      
end   % for frameCount

% Set the mouse pointer to normal again
set(gcf,'Pointer','arrow');

% Here is where all the plotting starts
if ptPostpro.neighbourplot_1

    % Generate the avg velocity plot (all cells)
    h_fig = figure('Name', imageName);

    % Draw a plot showing average velocity of all cells
    ymax = max (avgTrajFrame) + 1;
    plot (xAxis, avgTrajFrame); 
    title ('Avg Neighbour Trajectory Length');
    xlabel ('Frames');
    ylabel ('Avg Traj. Length (um)');
    if length (xAxis) > 1
       axis ([xAxis(1) xAxis(end) 0 ymax]);
    else
       axis ([xAxis(1) xAxis(1)+1 0 ymax]);
    end

    % Save the figures in fig, eps and tif format        
    hgsave (h_fig,[savePath filesep 'avgNeighbourTrajLength.fig']);
    print (h_fig, [savePath filesep 'avgNeighbourTrajLength.eps'],'-depsc2','-tiff');
    print (h_fig, [savePath filesep 'avgNeighbourTrajLength.tif'],'-dtiff');      
    
     % Save CSV files
    csvwrite ('avgNeighbourTrajLength.csv', [xAxis ; avgTrajFrame]);
end
