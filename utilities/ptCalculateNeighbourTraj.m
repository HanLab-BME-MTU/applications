function [neighTrajStats, xAxis] = ptCalculateNeighbourTraj (handles)
% ptCalculateNeighbourTraj plots neighbour traj. information gathered in MPM. 
%
% SYNOPSIS       ptCalculateNeighbourTraj (ptPostpro, MPM)
%
% INPUT          handles : a structure which contains the information from the GUI
%                MPM       : matrix containing the cell tracks
%                
% OUTPUT         neighTrajStats : struct with following fields:
%                    avgTrajFrame : vector with avg trajectory length changes
%
% DEPENDENCIES   ptCalculateNeighbourTraj  uses {nothing}
%                                  
%                ptCalculateNeighbourTraj is used by { PolyTrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Jul 04          Initial version
% Andre Kerstens        Aug 04          Renamed to ptCalculateNeighbourTraj
% Andre Kerstens        Sep 04          Complete rewrite of plotting functions

% Get the latest data from the handles
MPM = handles.allMPM;
cellProps = handles.allCellProps;
clusterProps = handles.allClusterProps;
frameProps = handles.allFrameProps;
validFrames = handles.allValidFrames;
jobData = handles.jobData;
guiData = handles.guiData;

% Get values from the gui (these are used for all jobs)
plotStartFrame = guiData.plotfirstimg;
plotEndFrame = guiData.plotlastimg;
maxDistance = guiData.maxneighbourdist;

% Determine the movie with the most frames
%[longestMPM, mpmLength] = ptMaxMPMLength (MPM);
%maxFrames = mpmLength / 2;

% Determine the movie with the most frames
%[shortestMPM, mpmLength] = ptMinMPMLength (MPM);
%minFrames = mpmLength / 2;
[shortestMovie, minFrames] = ptMinMovieLength (validFrames);

% Make sure we only process up to the shortest MPM else the averaging will
% not work correctly
if plotEndFrame > minFrames
    plotEndFrame = minFrames;
end

% Get start and end frames and increment value
%startFrame = jobData(1).firstimg;
%endFrame = jobData(shortestMPM).lastimg;
%increment = jobData(1).increment;
%numberOfFrames = ceil((plotEndFrame - plotStartFrame) / increment) + 1;
startFrame = jobData(shortestMovie).firstimg;
endFrame = jobData(shortestMovie).lastimg;
increment = jobData(shortestMovie).increment;
numberOfFrames = ceil((plotEndFrame - plotStartFrame) / increment) + 1;

% Get pixellength and frame interval
frameInterval = round (jobData(1).timeperframe / 60);    % In minutes
pixelLength = jobData(1).mmpixel;

% Over how many frames to we calculate trajectories
% eg 5 frames means 4 trajectories
nrOfTrajectories = guiData.nrtrajectories - 1;

% Initialize the trajectory matrix
avgTrajFrame = zeros (1, numberOfFrames-nrOfTrajectories-1);

% Initialize properties counter depending on radiobutton value
alwaysCountFrom1 = get (handles.GUI_alwayscount1_cb, 'Value');
if ~alwaysCountFrom1
   MPMCount = ceil ((plotStartFrame - startFrame) / increment) + 1;
else
   MPMCount = ceil ((plotStartFrame - 1) / increment) + 1;
end

% Initialize X-axis vector
xAxis = zeros (1, numberOfFrames-nrOfTrajectories-1);
iCount = 0;

% Go through every frame of the set.
for frameCount = plotStartFrame : increment : plotEndFrame
    
   % Increase MPM counter
   MPMCount = MPMCount + 1;
   
   % Initialize average sum counter
   averageSum = 0;

   % There should be enough frames left at the end to calc. traj.
   if (MPMCount + nrOfTrajectories) < numberOfFrames
   
      % Increase counter
      %iCount = iCount + 1;
       
      % Initialize counter for cells with nearest-enough neighbours
      cellsWithNeighbours = 0;

      % Clear variables
      clear cells triangleIndex cellEntries neighbours;
      clear distance mpmIndex avgTrajNeighbours;
      clear trajFrame;
      
      for jobCount = 1 : length(MPM)   
          
          % Find the index where this frame can be found
          frameIndx = find(validFrames{jobCount}(1,:) == MPMCount);

          if isempty(frameIndx)
              % Frame was bad and cannot be found in MPM
              avgTrajNeighbours{jobCount} = [];
          else
              
              % Increase counter
              if jobCount == 1
                 iCount = iCount + 1;
              end
      
              % Get the cell list for this frame and throw the zero entries out
              cells{jobCount} = MPM{jobCount}(:, 2*frameIndx-1 : 2*frameIndx);
              cells{jobCount} = cells{jobCount}(find (cells{jobCount}(:,1) ~= 0 & cells{jobCount}(:,2) ~= 0),:);

              % Triangulate all the cells with their neighbours
              triangleIndex{jobCount} = delaunay (cells{jobCount}(:,1), cells{jobCount}(:,2));

              % Initialize avgTraj
              avgTrajNeighbours{jobCount} = zeros (1, size(cells{jobCount},1));


              % Find neighbours for all these cells and do some trajectory calculations
              for cCount = 1 : size(cells{jobCount},1)

                 % Find the entries in triangleIndex for cell 'iCount'
                 cellEntries{jobCount} = triangleIndex{jobCount}(find (triangleIndex{jobCount}(:,1)==cCount | ...
                                                                       triangleIndex{jobCount}(:,2)==cCount | ...
                                                                       triangleIndex{jobCount}(:,3)==cCount),:);

                 % Extract the neighbours (incl duplicates)
                 nCount = 1;
                 for hCount = 1 : size (cellEntries{jobCount}, 1)
                    if cellEntries{jobCount}(hCount,1) == cCount
                       neighbours{jobCount}(nCount,1) = cellEntries{jobCount}(hCount,2);
                       nCount = nCount + 1;
                       neighbours{jobCount}(nCount,1) = cellEntries{jobCount}(hCount,3);
                       nCount = nCount + 1;
                    end
                    if cellEntries{jobCount}(hCount,2) == cCount
                       neighbours{jobCount}(nCount,1) = cellEntries{jobCount}(hCount,1);
                       nCount = nCount + 1;
                       neighbours{jobCount}(nCount,1) = cellEntries{jobCount}(hCount,3);
                       nCount = nCount + 1;
                    end
                    if cellEntries{jobCount}(hCount,3) == cCount
                       neighbours{jobCount}(nCount,1) = cellEntries{jobCount}(hCount,1);
                       nCount = nCount + 1;
                       neighbours{jobCount}(nCount,1) = cellEntries{jobCount}(hCount,2);
                       nCount = nCount + 1;
                    end
                 end
                 % Take only the unique neighbours
                 neighbours{jobCount} = unique (neighbours{jobCount});

                 % Calculate the distance to all of the neighbours
                 distance{jobCount} = [];
                 for dCount = 1 : length (neighbours{jobCount})   
                    distance{jobCount}(dCount,1) = sqrt ((cells{jobCount}(cCount, 1) - ...
                                                          cells{jobCount}(neighbours{jobCount}(dCount,1), 1))^2 + ...
                                                         (cells{jobCount}(cCount, 2) - ...
                                                          cells{jobCount}(neighbours{jobCount}(dCount,1), 2))^2);  
                 end

                 % Only keep the neighbours that are close enough
                 neighbours{jobCount} = neighbours{jobCount}(find (distance{jobCount} < maxDistance));

                 % In case the cell has any neighbours do some work
                 if ~isempty (neighbours{jobCount})

                    cellsWithNeighbours = cellsWithNeighbours + 1;

                    % Initialize avgTraj
                    avgTraj = zeros (1, length (neighbours{jobCount}));

                    % For every neighbour we are going to calculate avg trajectory
                    % length for the specified amount of frames
                    for nCount = 1 : length (neighbours{jobCount})

                       % Fetch the neighbour from cells
                       curNeighbour = cells{jobCount}(neighbours{jobCount}(nCount,1),:);

                       % Find the current neighbour in MPM
                       mpmIndex = find (MPM{jobCount}(:,2*frameIndx-1) == curNeighbour(:,1) & ...
                                        MPM{jobCount}(:,2*frameIndx) == curNeighbour(:,2));

                       % If multiple are found (for whatever reason), pick
                       % the first one
                       if length(mpmIndex) > 1
                           mpmIndex = mpmIndex(1);
                       end
                                        
                       % Calculate trajectories for the specified amount of frames for this neighbour
                       traj = [];
                       for tCount = 0 : nrOfTrajectories
                          traj(tCount+1,1) = sqrt((MPM{jobCount}(mpmIndex, 2*frameIndx-1+(tCount*2)) - ...
                                                   MPM{jobCount}(mpmIndex, 2*frameIndx-1+(tCount*2+2)))^2 + ...
                                                  (MPM{jobCount}(mpmIndex, 2*frameIndx+(tCount*2)) - ...
                                                   MPM{jobCount}(mpmIndex, 2*frameIndx+(tCount*2+2)))^2);
                          % Convert into micrometer per minute
                          traj(tCount+1) = (traj(tCount+1) * pixelLength) / (validFrames{jobCount}(2,frameIndx)/60);
                       end

                       % Average the trajectory
                       avgTraj(nCount) = (sum(traj) / (nrOfTrajectories+1));

                    end  % for nCount = 1 : length (neighbours{jobCount})

                    % Average the trajectory for all neighbours
                    avgTrajNeighbours{jobCount}(cCount) = sum(avgTraj) / length (neighbours{jobCount}); 

                 end  % if ~isempty (neighbours{jobCount})
              end  % for cCount = 1 : size(cells{jobCount},1)
              
              % Increase counter used later to calculate average
              averageSum = averageSum + 1;
              
          end  % if isempty(frameIndx)
      end  % jobCount = 1 : length(MPM)   

      if averageSum > 0   % Only go on if we have at least 1 good frame in the joblist
      
          % Store the frame number for display on the x-axis
          if ~alwaysCountFrom1
              xAxis(iCount) = frameCount;
          else
              xAxis(iCount) = iCount;
          end
          
          % Cat the avgTrajNeighbours together
          catAvgTrajNeighbours = cat (2, avgTrajNeighbours{:});

          % Sum and average the trajectory frames
          avgTrajFrame(iCount) = sum (catAvgTrajNeighbours) / cellsWithNeighbours;
          
      end  % if averageSum > 0
   end   % if (MPMCount + nrOfTrajectories) < numberOfFrame
end   % for frameCount

% Prepare output data
neighTrajStats.avgTrajFrame = avgTrajFrame(1:iCount);

% Make sure the x-axis has the correct length
xAxis = xAxis(1:iCount);