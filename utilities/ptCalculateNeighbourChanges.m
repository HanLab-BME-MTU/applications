function ptCalculateNeighbourChanges (ptPostpro, MPM)
% ptCalculateNeighbourChanges plots neighbour change info from MPM. 
%
% SYNOPSIS       ptCalculateNeighbourChanges (ptPostpro, MPM)
%
% INPUT          ptPostpro : a structure which contains the information
%                            from the GUI
%                MPM       : matrix containing the cell tracks
%                
% OUTPUT         None (plots are directly shown on the screen) 
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
maxDistance = ptPostpro.maxdistance;
maxDistToNeighbour = maxDistance;

% Initialize the avg neighbour change vector
avgNbChange = zeros (1, numberOfFrames);

% Initialize MPM counter (start with 0 if plotStartFrame=startFrame)
MPMCount = ceil ((plotStartFrame - startFrame) / increment);

% Initialize X-axis vector
xAxis = zeros (1, numberOfFrames-1);
iCount = 0;

% These calculations can take a while so set the mouse pointer to busy
set(gcf,'Pointer','watch');

% Go through every frame of the set.
for frameCount = plotStartFrame : increment : plotEndFrame
    
   % Increase MPM counter
   MPMCount = MPMCount + 1;
   
   % Store the frame number for display on the x-axis
   iCount = iCount + 1;
   xAxis (iCount) = frameCount;
   
   % Get the cell list for this frame, throw the zero entries out for
   % the neighbour change calculations, but keep the original indexnrs
   cells = MPM (:, 2*MPMCount-1 : 2*MPMCount);
   [cellIndexRow, cellIndexCol] = find (cells);
   cellIndex = unique(cellIndexRow)';
      
   % Triangulate all the cells with their neighbours
   triangleIndex = delaunay (cells(:,1), cells(:,2));
      
   % Initialize counter for cells with nearest-enough neighbours
   cellsWithNeighbours = 0;
      
   % Initialize neighbour change counter
   nbChange = zeros(length(cellIndex),1);
   
   % Initialize counters
   count = 1;
   neighCount = 1;
   
   % Init neighbours vector
   clear neighbours;
   neighbours(1:length(cellIndex))=struct('cell',[],'neighbours',[]);
   
   % Find neighbours for all these cells and do neighbour change analysis
   for cCount = cellIndex
                 
      % Find the entries in triangleIndex for cell 'cCount'
      neighboursTemp = [];
      cellEntries = triangleIndex(find (triangleIndex(:,1)==cCount | triangleIndex(:,2)==cCount | ...
                                        triangleIndex(:,3)==cCount),:);
         
      % Continue if there are no cells in the triangualtion matrix                                  
      if isempty(cellEntries)
         continue;
      end
                                        
      % Extract the neighbours (incl duplicates)
      nCount = 1;
      for hCount = 1 : size (cellEntries, 1)
         if cellEntries(hCount,1) == cCount
            neighboursTemp(nCount,1) = cellEntries(hCount,2);
            nCount = nCount + 1;
            neighboursTemp(nCount,1) = cellEntries(hCount,3);
            nCount = nCount + 1;
         end
         if cellEntries(hCount,2) == cCount
            neighboursTemp(nCount,1) = cellEntries(hCount,1);
            nCount = nCount + 1;
            neighboursTemp(nCount,1) = cellEntries(hCount,3);
            nCount = nCount + 1;
         end
         if cellEntries(hCount,3) == cCount
            neighboursTemp(nCount,1) = cellEntries(hCount,1);
            nCount = nCount + 1;
            neighboursTemp(nCount,1) = cellEntries(hCount,2);
            nCount = nCount + 1;
         end
      end
      
      % Take only the unique neighbours
      neighboursTemp = unique (neighboursTemp);
      
      % Find the 0-entry triangle point (last entry in cells)
      unrealNeighbour = length (cells);
      unrealNeighIndex = find (neighboursTemp == unrealNeighbour);
      
      % If it exists, kick it out of the neighbour vector
      if ~isempty(unrealNeighIndex)
         neighboursTemp(unrealNeighIndex) = [];
      end
         
      % Calculate the distance to all of the neighbours
      distance = [];
      for dCount = 1 : length (neighboursTemp)   
         distance(dCount,1) = sqrt ((cells(cCount, 1) - cells(neighboursTemp(dCount,1), 1))^2 + ...
                                    (cells(cCount, 2) - cells(neighboursTemp(dCount,1), 2))^2);  
      end
        
      % Only keep the neighbours that are close enough
      neighbours(neighCount).cell = cCount;
      neighbours(neighCount).neighbours = neighboursTemp(find (distance < maxDistToNeighbour))';                                 

      % To compare neighbours we need to be at least at frame 2
      if frameCount > plotStartFrame
           
         % Get the cell numbers from the previous neighbours cells
         %prevCellNrs = [prevNeighbours.cell]';
         
         % Find the index for the current cell in prevCellNrs
         prevNeighbourInd = find(prevCellNrs == cCount);
         
         % In case the current cellnr can be found in the previous
         % neighbours we go ahead and compare
         if ~isempty(prevNeighbourInd)
             
            % Get the index for neighbours as well
            cellNrs = [neighbours.cell]';
            neighbourInd = find (cellNrs == cCount);
         
            % Both contain the cell: find out if neighbours are similar
            if length(neighbours(neighbourInd).neighbours) >= length(prevNeighbours(prevNeighbourInd).neighbours)    
               neighbourChange = ~ismember (neighbours(neighbourInd).neighbours, prevNeighbours(prevNeighbourInd).neighbours); 
            else
               neighbourChange = ~ismember (prevNeighbours(prevNeighbourInd).neighbours, neighbours(neighbourInd).neighbours); 
            end
            
            % Sum up the changes
            nbChange(count) = nbChange(count) + sum(neighbourChange);
            
            % Increase change counter
            count = count + 1;
         end  % ~isempty(find(prevCellNrs == cCount))
      end  % if frameCount > plotStartFrame 
      
      % Increase neighbour struct counter
      neighCount = neighCount + 1;
      
   end  % for cCount = cellIndex
   
   % Average the neighbour changes
   avgNbChange(iCount) = sum(nbChange); %/ length(nbChange);
      
   % Keep the neighbours of this frame for use in the next
   prevNeighbours = neighbours;
   
   % Store the cell numbers from the previous neighbours cells
   prevCellNrs = [prevNeighbours.cell]';
end   % for frameCount

% Set the mouse pointer to normal again
set(gcf,'Pointer','arrow');

if ptPostpro.neighbourplot_2

    % Generate the neighbour change plot (all cells)
    h_fig2 = figure('Name', imageName);

    % Draw a plot showing average velocity of all cells
    ymax = max (avgNbChange) + (0.1*max (avgNbChange));
    plot (xAxis, avgNbChange); 
    title ('Avg Neighbour Interaction Change');
    xlabel ('Frames');
    ylabel ('Avg Neighbour Change');
    if length (xAxis) > 1
       axis ([xAxis(1) xAxis(end) 0 ymax]);
    else
       axis ([xAxis(1) xAxis(1)+1 0 ymax]);
    end

    % Save the figures in fig, eps and tif format        
    hgsave (h_fig2,[savePath filesep 'avgNeighbourChange.fig']);
    print (h_fig2, [savePath filesep 'avgNeighbourChange.eps'],'-depsc2','-tiff');
    print (h_fig2, [savePath filesep 'avgNeighbourChange.tif'],'-dtiff');      
    
     % Save CSV files
    csvwrite ('avgNeighbourInteractionChange.csv', [xAxis ; avgNbChange]);
end

