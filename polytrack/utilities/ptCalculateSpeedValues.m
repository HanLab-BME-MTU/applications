function [avgVelocityStats, velocitySingleStats, velocityVarStats, velocityHistStats, xAxis] = ptCalculateSpeedValues (handles)
% ptPlotSpeedValues plots speed information gathered in MPM. 
%
% SYNOPSIS       [avgVelocityStats, velocitySingleStats, velocityVarStats, velHistStats, xAxis] = 
%                               ptCalculateSpeedValues (fileList, filesSelected)
%
% INPUT          handles : the gui handles struct
%                
% OUTPUT         avgVelocityStats : struct with following fields:
%                   avgVelocity : vector with avg velocity all cells
%                   avgVelocitySquared : vector with squared avg velocity all cells
%                   avgSingleVelocity : vector with avg single cell vel.
%                   avgClusteredVelocity : vector with avg clustered cell vel.
%                velocitySingleStats : struct with following fields:
%                   velAllCellsHigherThanAvgSingleCells : vector with
%                      amount of all cells with velocity higher than avg single cells
%                velocityVarStats : struct with following fields:
%                   varVelocity : vector with velocity variance all cells
%                   varSingleCellVelocity : vector with velocity variance single cells
%                   varClusteredCellVelocity : vector with velocity variance clustered cells
%                velHistStats : struct with following fields:
%                   velocityHist : velocity histogram values
%                   maxVelocity : vector with max velocities histogram
%                xAxis : vector with x-axis values
%
% DEPENDENCIES   ptCalculateSpeedValues  uses {nothing}
%                                  
%                ptCalculateSpeedValues is used by { PolyTrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Colin Glass           Feb 04          Initial release
% Andre Kerstens        Jun 04          Cleaned up source and renamed file
% Andre Kerstens        Jun 04          Changed ymax to 100% for percentage all cells
% Andre Kerstens        Aug 04          Fixed bug in velocity squared and made loop more 
%                                       effective (less calculations)
% Andre Kerstens        Sep 04          Complete rewrite of plot functionality

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
binSize = guiData.binsize;
multipleFrameVelocity = guiData.multframevelocity;

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
% startFrame = jobData(1).firstimg;
% endFrame = jobData(shortestMPM).lastimg;
% increment = jobData(1).increment;
% numberOfFrames = ceil((plotEndFrame - plotStartFrame) / increment) + 1;
startFrame = jobData(shortestMovie).firstimg;
endFrame = jobData(shortestMovie).lastimg;
increment = jobData(shortestMovie).increment;
numberOfFrames = ceil((plotEndFrame - plotStartFrame) / increment) + 1;

% Get pixellength and frame interval
% Note: we can get the framerate from the validFrames array
%frameInterval = round (jobData(1).timeperframe / 60);    % In minutes
pixelLength = jobData(1).mmpixel;

% Initialize the displacement and x-axis matrices
avgDisplacement = zeros (1, numberOfFrames-multipleFrameVelocity);
avgSingleDisplacement = zeros (1, numberOfFrames-multipleFrameVelocity);
avgClusteredDisplacement = zeros (1, numberOfFrames-multipleFrameVelocity);
avgVelocity = zeros (1, numberOfFrames-multipleFrameVelocity);
avgVelocitySquared = zeros (1, numberOfFrames-multipleFrameVelocity);
avgSingleVelocity = zeros (1, numberOfFrames-multipleFrameVelocity);
avgClusteredVelocity = zeros (1, numberOfFrames-multipleFrameVelocity);
velAllCellsHigherThanAvgSingleCells = zeros (1, numberOfFrames-multipleFrameVelocity);
velocityHist = zeros (binSize, numberOfFrames-multipleFrameVelocity);
maxVelocity = zeros (1, numberOfFrames-multipleFrameVelocity);
xAxis = zeros (1, numberOfFrames-multipleFrameVelocity);

% Initialize MPM counter
%MPMCount = ceil ((plotStartFrame - startFrame) / increment) + 1;

% Initialize properties counter depending on radiobutton value
alwaysCountFrom1 = get (handles.GUI_alwayscount1_cb, 'Value');
if ~alwaysCountFrom1
   MPMCount = ceil ((plotStartFrame - startFrame) / increment) + 1;
else
   MPMCount = ceil ((plotStartFrame - 1) / increment) + 1;
end

% Initialize x-axis and counter for avg vectors
xAxis = zeros (1, numberOfFrames-multipleFrameVelocity);
iCount = 0;

% Go through every frame of the set. Start at the second frame
% because only from there we can start calculating a displacement
for frameCount = plotStartFrame+increment : increment : plotEndFrame

   % Increase MPM counter
   MPMCount = MPMCount + 1;

   % Initialize average sum counter
   averageSum = 0;

   % If the velocity is calculated over multiple frames we should have
   % enough frames to do this
   if MPMCount > find(validFrames{shortestMovie}(1,:) == multipleFrameVelocity)

      % Increase counter
      %iCount = iCount + 1;

      for jobCount = 1 : length(MPM)

          % Find the index where this frame can be found
          frameIndx = find(validFrames{jobCount}(1,:) == MPMCount);

          if isempty(frameIndx) | frameIndx > length(cellProps{jobCount})
              % Frame was bad and cannot be found in MPM
              displacement{jobCount} = [];
              velocity{jobCount} = [];
              singleCellsDisplacement{jobCount} = [];
              singleCellVelocity{jobCount} = [];
              clusteredCellsDisplacement{jobCount} = [];
              clusteredCellVelocity{jobCount} = [];
          else
              
              % Increase counter
              if jobCount == 1
                 iCount = iCount + 1;
              end
          
              % Get the coordinates from the previous (dependent on the value of
              % multipleFrameVelocity) and the current frame
              %prevCoordinates{jobCount} = MPM{jobCount}(:, 2*(MPMCount-multipleFrameVelocity)-1 : 2*(MPMCount-multipleFrameVelocity));
              %curCoordinates{jobCount} = MPM{jobCount}(:, 2*MPMCount-1 : 2*MPMCount);
              %coordinates{jobCount} = [prevCoordinates{jobCount} curCoordinates{jobCount}];
              prevCoordinates{jobCount} = MPM{jobCount}(:, 2*(frameIndx-multipleFrameVelocity)-1 : 2*(frameIndx-multipleFrameVelocity));
              curCoordinates{jobCount} = MPM{jobCount}(:, 2*frameIndx-1 : 2*frameIndx);
              coordinates{jobCount} = [prevCoordinates{jobCount} curCoordinates{jobCount}];

              % If coordinates in the current frame have started a track or the
              % ones in the previous frame ended a track, we cannot
              % calculate a displacement so we should throw these out
              coordinates{jobCount}(find ((coordinates{jobCount}(:,1) == 0 & coordinates{jobCount}(:,2) == 0) | ...
                                          (coordinates{jobCount}(:,3) == 0 & coordinates{jobCount}(:,4) == 0)),:) = [];

              % From the cluster and cell properties we can find out which cells
              % are single and which belong to a cluster. First initialize two
              % matrices for this
              singleCellIndex{jobCount} = zeros (size (coordinates{jobCount},1), 1);
              clusteredCellIndex{jobCount} = zeros (size (coordinates{jobCount},1), 1);

              % Initialize single and clustered cell coordinates counters
              singleCellCount = 0;
              clusteredCellCount = 0;

              % Find the single cell and clustered coordinates and store these
              for jCount = 1 : size (coordinates{jobCount},1) 

                 % Find the index of the cell in CellProps
                 cellIndex = find (cellProps{jobCount}(:,1,frameIndx) == coordinates{jobCount}(jCount,3) & ...
                                   cellProps{jobCount}(:,2,frameIndx) == coordinates{jobCount}(jCount,4));

                 % In case the cell wasn't found in cellProps (maybe it was on the background?)
                 % we can skip over this part
                 if ~isempty (cellIndex)

                    % Check whether cellIndex contains double coordinate
                    % pairs (shouldn't be)
                    if length(cellIndex) > 1
                        cellIndex = cellIndex(1);   % Take the first one
                    end
                     
                    % Find this cluster in clusterProps
                    clusterIndex = find (clusterProps{jobCount}(:,1,frameIndx) == cellProps{jobCount}(cellIndex,3,frameIndx));

                    % Test whether this cell is in a cluster or not
                    if clusterProps{jobCount}(clusterIndex,2,frameIndx) == 1  % Single cell found

                       % Increase the counter for the single cell matrix
                       singleCellCount = singleCellCount + 1; 

                       % Store the coordinates
                       singleCellIndex{jobCount}(singleCellCount,:) = jCount;
                    else    % The cell was in a cluster

                       % Increase the counter for the clustered cell matrix
                       clusteredCellCount = clusteredCellCount + 1; 

                       % Store the coordinates
                       clusteredCellIndex{jobCount}(clusteredCellCount,:) = jCount;
                    end
                 end  % if ~isempty
              end  % for jCount = 1 : size (coordinates,1)

              % Remove all the zero entries in both matrices
              singleCellIndex{jobCount}(find (singleCellIndex{jobCount}(:,1) == 0), :) = [];
              clusteredCellIndex{jobCount}(find (clusteredCellIndex{jobCount}(:,1) == 0), :) = [];  

              % Calculate the displacement for all cells
              displacement{jobCount} = sqrt ((coordinates{jobCount}(:,1) - coordinates{jobCount}(:,3)).^2 + ...
                                             (coordinates{jobCount}(:,2) - coordinates{jobCount}(:,4)).^2);

              % Calculate true velocity in um/min for all cells; use the
              % validFrames array for this
              velocity{jobCount} = (displacement{jobCount} * pixelLength) / ...
                                   (multipleFrameVelocity * (validFrames{jobCount}(2,frameIndx)/60));

              % Calculate the displacement for single cells
              singleCellsDisplacement{jobCount} = displacement{jobCount}(singleCellIndex{jobCount},:);

              % Calculate true velocity in um/min for all cells
              singleCellVelocity{jobCount} = velocity{jobCount}(singleCellIndex{jobCount},:);

              % Calculate the displacement for clustered cells
              clusteredCellsDisplacement{jobCount} = displacement{jobCount}(clusteredCellIndex{jobCount},:);

              % Calculate true velocity in um/min for all cells
              clusteredCellVelocity{jobCount} = velocity{jobCount}(clusteredCellIndex{jobCount},:);
              
              % Increase counter used later to calculate average
              averageSum = averageSum + 1;

          end  % if isempty(frameIndx)
      end  % for jobCount = 1 : length(MPM)
        
      if averageSum > 0   % Only go on if we have at least 1 good frame in the joblist

          % Store x-Axis value
          if ~alwaysCountFrom1
              xAxis(iCount) = MPMCount;
          else
              xAxis(iCount) = iCount;
          end

          % Cat all the matrices that we found together
          allDisplacement = cat(1,displacement{:});
          allVelocity = cat(1, velocity{:});
          allSingleCellsDisplacement = cat(1, singleCellsDisplacement{:});
          allSingleCellVelocity = cat(1, singleCellVelocity{:});
          allClusteredCellsDisplacement = cat(1, clusteredCellsDisplacement{:});
          allClusteredCellVelocity = cat(1, clusteredCellVelocity{:});

          % Calculate histogram for velocity all cells
          histVelAllCells{iCount} = hist (allVelocity, binSize);
          %filename = [SaveDir filesep 'histVelAllCells_' num2str(iCount) '.mat'];
          %save (filename, 'histVelAllCells');

          % Save the velocity hist for single cells
          histVelSingleCells{iCount} = hist (allSingleCellVelocity, binSize);
          %filename = [SaveDir filesep 'histVelSingleCells_' num2str(iCount) '.mat'];
          %save (filename, 'histVelSingleCells');

          % Save the velocity hist for clustered cells
          histVelClusteredCells{iCount} = hist (allClusteredCellVelocity, binSize);
          %filename = [SaveDir filesep 'histVelClusteredCells_' num2str(iCount) '.mat'];
          %save (filename, 'histVelClusteredCells');

          % Store max velocity of the frame
          if ~isempty(allVelocity)
             maxVelocity(iCount) = max (allVelocity);
          else
             maxVelocity(iCount) = 0;
          end

          % Store max single cell velocity of the frame
          if ~isempty(allSingleCellVelocity)
             maxSingleCellVelocity(iCount) = max (allSingleCellVelocity);
          else
             maxSingleCellVelocity(iCount) = 0;
          end

          % Store max clustered cell velocity of the frame
          if ~isempty(allClusteredCellVelocity)
             maxClusteredCellVelocity(iCount) = max (allClusteredCellVelocity);
          else
             maxClusteredCellVelocity(iCount) = 0;
          end      

          % Calculate the variance of the velocity for all cells
          if ~isempty (allVelocity)
             varVelocity(iCount) = var (allVelocity);
          else
             varVelocity(iCount) = 0;
          end

          % Calculate the variance of the velocity for single cells
          if ~isempty (allSingleCellVelocity)
             varSingleCellVelocity(iCount) = var (allSingleCellVelocity);
          else
             varSingleCellVelocity(iCount) = 0;
          end

          % Calculate the variance of the velocity for clustered cells
          if ~isempty (allClusteredCellVelocity)
             varClusteredCellVelocity(iCount) = var (allClusteredCellVelocity);
          else
             varClusteredCellVelocity(iCount) = 0;
          end

          % From this the average displacement for all cells can be calculated
          if length (allDisplacement) > 0
             avgDisplacement(iCount) = sum (allDisplacement) / length (allDisplacement);
          else
             avgDisplacement(iCount) = 0;
          end

          % From this the average displacement for single cells can be calculated
          if length (allSingleCellsDisplacement) > 0
             avgSingleDisplacement(iCount) = sum (allSingleCellsDisplacement) / length (allSingleCellsDisplacement);
          else
             avgSingleDisplacement(iCount) = 0;
          end

          % From this the average displacement for clustered cells can be calculated
          if length (allClusteredCellsDisplacement) > 0
             avgClusteredDisplacement(iCount) = sum (allClusteredCellsDisplacement) / length (allClusteredCellsDisplacement);
          else
             avgClusteredDisplacement(iCount) = 0;
          end

          % From this the average velocity for all cells can be calculated
          if length (allVelocity) > 0
             avgVelocity(iCount) = sum (allVelocity) / length (allVelocity);
             avgVelocitySquared(iCount) = avgVelocity(iCount)^2;
          else
             avgVelocity(iCount) = 0;
             avgVelocitySquared(iCount) = 0;
          end

          % From this the average velocity for single cells can be calculated
          if length (allSingleCellsDisplacement) > 0
             avgSingleVelocity (iCount) = sum (allSingleCellVelocity) / length (allSingleCellVelocity);
          else
             avgSingleVelocity (iCount) = 0;
          end

          % From this the average velocity for clustered cells can be calculated
          if length (allClusteredCellsDisplacement) > 0
             avgClusteredVelocity(iCount) = sum (allClusteredCellVelocity) / length (allClusteredCellVelocity);
          else
             avgClusteredVelocity(iCount) = 0;
          end

          % Calculate how many cells (of all cells) have a velocity higher or
          % equal than the average single cell velocity (in percent)
          if length (allVelocity) ~= 0
             velAllCellsHigherThanAvgSingleCells(iCount) = (length (find (allVelocity >= avgSingleVelocity(iCount))) / ...
                                                            length (allVelocity)) * 100.0;
          else
             velAllCellsHigherThanAvgSingleCells(iCount) = 0;
          end     

          % Generate a velocity histogram
          allVelocityHist (:,iCount) = hist (allVelocity, binSize);

       end   % if averageSum  
   end   % if MPMCount > multipleFrameVelocity
end   % for frameCount

% Store everything in a struct for easier handling outside of this function
avgVelocityStats.avgVelocity = avgVelocity(1:iCount);
avgVelocityStats.avgVelocitySquared = avgVelocitySquared(1:iCount);
avgVelocityStats.avgSingleVelocity = avgSingleVelocity(1:iCount);
avgVelocityStats.avgClusteredVelocity = avgClusteredVelocity(1:iCount);
avgVelocityStats.avgSingleDisplacement = avgSingleDisplacement(1:iCount);

velocitySingleStats.velAllCellsHigherThanAvgSingleCells = velAllCellsHigherThanAvgSingleCells(1:iCount);

velocityVarStats.varVelocity = varVelocity(1:iCount);
velocityVarStats.varSingleCellVelocity = varSingleCellVelocity(1:iCount);
velocityVarStats.varClusteredCellVelocity = varClusteredCellVelocity(1:iCount);

velocityHistStats.velocityHist = allVelocityHist;
velocityHistStats.histVelAllCells = histVelAllCells;
velocityHistStats.histVelSingleCells = histVelSingleCells;
velocityHistStats.histVelClusteredCells = histVelClusteredCells;
velocityHistStats.maxVelocity = maxVelocity;
velocityHistStats.maxSingleCellVelocity = maxSingleCellVelocity;
velocityHistStats.maxClusteredCellVelocity = maxClusteredCellVelocity;
velocityHistStats.binSize = guiData.binsize;

% Make sure the x-axis has the correct length
xAxis = xAxis(1:iCount);