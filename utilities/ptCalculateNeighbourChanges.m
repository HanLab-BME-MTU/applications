function [neighChangeStats, xAxis] = ptCalculateNeighbourChanges (handles)
% ptCalculateNeighbourChanges plots neighbour change info from MPM. 
%
% SYNOPSIS       ptCalculateNeighbourChanges (ptPostpro, MPM)
%
% INPUT          handles : a structure which contains the information from the GUI
%                
% OUTPUT         neighChangeStats : struct with following fields:
%                     : vector with avg trajectory length changes
%
% DEPENDENCIES   ptCalculateNeighbourChanges  uses {nothing}
%                                  
%                ptCalculateNeighbourhoodChanges is used by { PolyTrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Jul 04          Initial version of neighbour change calculation
% Andre Kerstens        Aug 04          Changed way how changes are counted; more accurate with new method
% Andre Kerstens        Aug 04          Changed neighbour datastore from cell array to struct (cell2mat 
%                                       function is way to slow!)
% Andre Kerstens        Aug 04          Fixed bug where csv where saved in wrong dir
% Andre Kerstens        Sep 04          Complete rewrite of plotting functions

% Get the latest data from the handles
MPM = handles.allMPM;
cellProps = handles.allCellProps;
clusterProps = handles.allClusterProps;
frameProps = handles.allFrameProps;
jobData = handles.jobData;
guiData = handles.guiData;

% Get values from the gui (these are used for all jobs)
plotStartFrame = guiData.plotfirstimg;
plotEndFrame = guiData.plotlastimg;
binSize = guiData.binsize;
multipleFrameVelocity = guiData.multFrameVelocity;
maxDistance = guiData.maxneighbourdist;

% Determine the movie with the most frames
[longestMPM, mpmLength] = ptMaxMPMLength (MPM);
maxFrames = mpmLength / 2;

% Get start and end frames and increment value
startFrame = jobData(1).firstimg;
endFrame = jobData(longestMPM).lastimg;
increment = jobData(1).increment;
numberOfFrames = ceil((plotEndFrame - plotStartFrame) / increment) + 1;

% Get pixellength and frame interval
frameInterval = round (jobData(1).timeperframe / 60);    % In minutes
pixelLength = jobData(1).mmpixel;

% Initialize the avg neighbour change vector
allAvgNbChange = zeros (1, numberOfFrames);

% Initialize MPM counter (start with 0 if plotStartFrame=startFrame)
MPMCount = ceil ((plotStartFrame - startFrame) / increment);

% Initialize X-axis vector
xAxis = zeros (1, numberOfFrames-1);
iCount = 0;

% Inform the user
%fprintf (1, '\nProcessing frame: ');

% Go through every frame of the set.
for frameCount = plotStartFrame : increment : plotEndFrame
    
   % Let the user know where we are
   %fprintf (1, '%d ', frameCount);
    
   % Increase MPM counter
   MPMCount = MPMCount + 1;
   
   % Store the frame number for display on the x-axis
   iCount = iCount + 1;
   xAxis (iCount) = frameCount;
   
   % Clear the some temp vars
   clear cells cellIndex neighbours neighboursTemp;
   clear nbChange avgNbChange;
   clear cellEntries cellNrs neighbourInd;
   clear cellIndexRow cellIndexCol unrealNeighIndex;
   clear prevNeighbourInd neighbourInd cellNrs;
   clear diff diffIndx;
   
   for jobCount = 1 : length(MPM) 
   
       % Get the cell list for this frame, throw the zero entries out for
       % the neighbour change calculations, but keep the original indexnrs
       cells{jobCount} = MPM{jobCount}(:, 2*MPMCount-1 : 2*MPMCount);
       [cellIndexRow{jobCount}, cellIndexCol{jobCount}] = find (cells{jobCount});
       cellIndex{jobCount} = unique(cellIndexRow{jobCount})';

       % Triangulate all the cells with their neighbours
       triangleIndex{jobCount} = delaunay (cells{jobCount}(:,1), cells{jobCount}(:,2));

       % Initialize counter for cells with nearest-enough neighbours
       cellsWithNeighbours = 0;

       % Initialize neighbour change counter
       nbChange{jobCount} = zeros(length(cellIndex{jobCount}),1);
       neighbourChangeCount = 0;

       % Initialize counters
       count = 1;
       neighCount = 1;

       % Init neighbours vector
       %clear neighbours;
       neighbours{jobCount}(1:length(cellIndex{jobCount})) = struct('cell',[],'neighbours',[]);

       % Find neighbours for all these cells and do neighbour change analysis
       for cCount = cellIndex{jobCount}

          % Find the entries in triangleIndex for cell 'cCount'
          neighboursTemp{jobCount} = [];
          cellEntries{jobCount} = triangleIndex{jobCount}(find (triangleIndex{jobCount}(:,1) == cCount | ...
                                                                triangleIndex{jobCount}(:,2) == cCount | ...
                                                                triangleIndex{jobCount}(:,3) == cCount),:);

          % Continue if there are no cells in the triangualtion matrix                                  
          if isempty(cellEntries{jobCount})
             continue;
          end

          % Extract the neighbours (incl duplicates)
          nCount = 1;
          for hCount = 1 : size (cellEntries{jobCount}, 1)
             if cellEntries{jobCount}(hCount,1) == cCount
                neighboursTemp{jobCount}(nCount,1) = cellEntries{jobCount}(hCount,2);
                nCount = nCount + 1;
                neighboursTemp{jobCount}(nCount,1) = cellEntries{jobCount}(hCount,3);
                nCount = nCount + 1;
             end
             if cellEntries{jobCount}(hCount,2) == cCount
                neighboursTemp{jobCount}(nCount,1) = cellEntries{jobCount}(hCount,1);
                nCount = nCount + 1;
                neighboursTemp{jobCount}(nCount,1) = cellEntries{jobCount}(hCount,3);
                nCount = nCount + 1;
             end
             if cellEntries{jobCount}(hCount,3) == cCount
                neighboursTemp{jobCount}(nCount,1) = cellEntries{jobCount}(hCount,1);
                nCount = nCount + 1;
                neighboursTemp{jobCount}(nCount,1) = cellEntries{jobCount}(hCount,2);
                nCount = nCount + 1;
             end
          end  % hCount = 1 : size (cellEntries, 1)

          % Take only the unique neighbours
          neighboursTemp{jobCount} = unique (neighboursTemp{jobCount});

          % Find the 0-entry triangle point (last entry in cells)
          unrealNeighbour = length (cells{jobCount});
          unrealNeighIndex{jobCount} = find (neighboursTemp{jobCount} == unrealNeighbour);

          % If it exists, kick it out of the neighbour vector
          if ~isempty(unrealNeighIndex)
             neighboursTemp{jobCount}(unrealNeighIndex{jobCount}) = [];
          end

          % Calculate the distance to all of the neighbours
          distance = [];
          for dCount = 1 : length (neighboursTemp{jobCount})   
             distance(dCount,1) = sqrt ((cells{jobCount}(cCount, 1) - ...
                                         cells{jobCount}(neighboursTemp{jobCount}(dCount,1), 1))^2 + ...
                                        (cells{jobCount}(cCount, 2) - ...
                                         cells{jobCount}(neighboursTemp{jobCount}(dCount,1), 2))^2);  
          end

          % Only keep the neighbours that are close enough
          neighbours{jobCount}(neighCount).cell = cCount;
          neighbours{jobCount}(neighCount).neighbours = neighboursTemp{jobCount}(find (distance < maxDistance))';                                 

          % To compare neighbours we need to be at least at frame 2
          if frameCount > plotStartFrame

             % Find the index for the current cell in prevCellNrs
             %prevNeighbourInd = find([prevNeighbours{jobCount}.cell]'prevCellNrs == cCount);
             prevNeighbourInd{jobCount} = find([prevNeighbours{jobCount}.cell]' == cCount);

             % In case the current cellnr can be found in the previous
             % neighbours we go ahead and compare
             if ~isempty(prevNeighbourInd{jobCount})

                % Get the index for neighbours as well
                cellNrs{jobCount} = [neighbours{jobCount}.cell]';
                neighbourInd{jobCount} = find (cellNrs{jobCount} == cCount);

                % Both contain the cell: find out if neighbours are similar
                if length(neighbours{jobCount}(neighbourInd{jobCount}).neighbours) >= ...
                       length(prevNeighbours{jobCount}(prevNeighbourInd{jobCount}).neighbours)    
                   % If a new cell comes into the image, we don't count it as a new neighbour
                   [diff, diffIndx] = setdiff(neighbours{jobCount}(neighbourInd{jobCount}).neighbours, ...
                                              prevNeighbours{jobCount}(prevNeighbourInd{jobCount}).neighbours);
                   for diffCount = 1 : length (diff)
                      if ismember(diff(diffCount), [prevNeighbours{jobCount}.cell]')
                      %if ismember(diff(diffCount), prevCellNrs)
                         % It is a real neighbour change, so count it
                         neighbourChangeCount = neighbourChangeCount + 1;
                      end
                   end
                else  % neighbours is shorter, which means we lost a cell(s)
                   % If the lost cell is due to bad tracking, we don't count it as a lost neighbour
                   [diff, diffIndx] = setdiff(prevNeighbours{jobCount}(prevNeighbourInd{jobCount}).neighbours, ...
                                              neighbours{jobCount}(neighbourInd{jobCount}).neighbours);
                   for diffCount = 1 : length (diff)
                      if ismember(diff(diffCount), cellNrs{jobCount})
                         % It is a real neighbour change, so count it
                         neighbourChangeCount = neighbourChangeCount + 1;
                      end
                   end
                end

                % Sum up the changes
                nbChange{jobCount}(count) = nbChange{jobCount}(count) + neighbourChangeCount;

                % Increase change counter
                count = count + 1;

                % Reset neighbourChangeCount
                neighbourChangeCount = 0;

             end  % ~isempty(find(prevCellNrs == cCount))
          end  % if frameCount > plotStartFrame 

          % Increase neighbour struct counter
          neighCount = neighCount + 1;

       end  % for cCount = cellIndex{jobCount}
       
       %avgNbChange{jobCount} = sum(nbChange{jobCount}) / length (cellIndex{jobCount});
       avgNbChange{jobCount} = sum(nbChange{jobCount});
       
   end  % for jobCount = 1 : length(MPM)
   
   % Cat the relevant vectors together
   catAvgNbChange = cat (1, avgNbChange{:})';
   catCellIndex = cat (2, cellIndex{:});
   
   % Average the neighbour changes and normalize over the number of cells
   allAvgNbChange(iCount) = sum(catAvgNbChange) / length (catCellIndex);
      
   % Keep the neighbours of this frame for use in the next
   prevNeighbours = neighbours;
end   % for frameCount

% Inform the user
%fprintf (1, '\n');

% Prepare output data
neighChangeStats.avgNbChange = allAvgNbChange;
