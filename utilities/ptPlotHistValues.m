function ptPlotHistValues (ptPostpro)
% ptPlotHistValues plots information gathered in cellProps and
% clusterProps regarding cell-cell distances 
%
% SYNOPSIS       ptPlotHistValues (ptPostpro)
%
% INPUT          ptPostpro : a structure which contains the information
%                            from the GUI (see below)
%                
% OUTPUT         None (plots are directly shown on the screen) 
%
% DEPENDENCIES   ptPlotHistValues.m  uses { nothing }
%                                  
%                ptPlotHistValues.m is used by { PolyTrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Jul 04          First version of ptPlotHistValues

% First assign all the postpro fields to a meaningfull variable
startFrame = ptPostpro.firstimg;
endFrame = ptPostpro.lastimg;
increment = ptPostpro.increment;
plotStartFrame = ptPostpro.plotfirstimg;
plotEndFrame = ptPostpro.plotlastimg;
numberOfFrames = ceil((plotEndFrame - plotStartFrame) / increment) + 1;
savePath = ptPostpro.saveallpath;
jobPath = ptPostpro.jobpath;
imageName = ptPostpro.imagename;
rowSize = ptPostpro.rowsize;
colSize = ptPostpro.colsize;
pixelLength = ptPostpro.mmpixel;
binSize = ptPostpro.binsize;
maxDistance = ptPostpro.maxdistance;

% Calculate the total area of the frame
totalAreaFrame = rowSize * colSize;

% Get radio button values
cellCellDistPlot = ptPostpro.cellcelldistplot;

% Get cell and cluster properties
cellProps = ptPostpro.cellProps;
clusterProps = ptPostpro.clusterProps;

% Initialize properties counter
propCount = ceil ((plotStartFrame - startFrame) / increment);

% Initialize average distance 
averageDist = [];

% Initialize X-axis vector and counter
xAxis = zeros (1, numberOfFrames);
iCount = 0;

% Calculate a number of statistics for every frame
for frameCount = plotStartFrame : increment : plotEndFrame
   
   % Update the properties counter
   propCount = propCount + 1;
    
   % Update the x-axis vector and counter
   iCount = iCount + 1;
   xAxis (iCount) = frameCount;
    
   % Remove the zero rows from cellProps 
   [notZeroEntryRows, notZeroEntryCols] = find (cellProps (:,:,propCount));
   notZeroEntryRows = unique (notZeroEntryRows);
   cells = cellProps (notZeroEntryRows,:,propCount);
   
   % Remove the zero rows from clusterProps 
   [notZeroEntryRows, notZeroEntryCols] = find (clusterProps (:,:,propCount));
   notZeroEntryRows = unique (notZeroEntryRows);
   clusters = clusterProps (notZeroEntryRows,:,propCount);
   
   % Calculate all the distances between cells by doing a Delaunay
   % triangulation
   triangles = delaunay (cells(:,1), cells(:,2));

   % Close the loop by copying the first coordinate set to a new column
   triangles(:,4) = triangles(:,1);

   % Match pairs of points, describing one line (no longer three points, describing a
   % triangle) and throw out the lines (between two points) that occur twice 
   uniqTriang = cat (1, unique(triangles(:,1:2), 'rows'), unique (triangles(:,2:3), 'rows'), ...
                        unique( triangles(:,3:4), 'rows'));
             
   % Prepare storage for the distances
   distances = zeros (length (uniqTriang), 1);

   % Calculate the distances
   for jCount = 1 : length (uniqTriang)   
      distances(jCount) = sqrt ((cells(uniqTriang(jCount,1), 1) - cells(uniqTriang(jCount,2), 1))^2 + ...
                                (cells(uniqTriang(jCount,1), 2) - cells(uniqTriang(jCount,2), 2))^2);  
   end

   % Calculate the average distance in micrometer
   averageDist(iCount) = (sum (distances) / length (distances)) * pixelLength;
   
   % Calculate a distance histogram
   %distanceHist (:,iCount) = hist (distances, binSize);
   
   
   if frameCount == plotStartFrame
      for kCount = 1 : size(cells,1)
         distAllCells = sqrt ((cells(kCount,1) - cells(:,1)).^2 + ...
                                       (cells(kCount,2) - cells(:,2)).^2);
   
         % Kick out all the ones that are > maxDistance
         neighbours = distAllCells (find (distAllCells < (3*maxDistance)));
         
         % Sort the neighbours vector
         sortedNeighbours = sort(neighbours);
         
         % Kick out the zero entry (the one distanced with itself)
         sortedNeighbours(find (sortedNeighbours == 0)) = [];
         
         % Take the 3 nearest neighbours, but leave out the zero entry
         if length (sortedNeighbours) >= 3
            nearNeighbours (kCount,:) = sortedNeighbours(1:3);
         else
            nearNeighbours(kCount,:) = [0 , 0 , 0];
            if length (sortedNeighbours) > 0
               nearNeighbours (kCount,1: length(sortedNeighbours)) = sortedNeighbours;
            end
         end
      end
   elseif frameCount > plotStartFrame
      % Dummy
   end
   
end 

if ptPostpro.cellcelldistplot_1

    % Generate the figure and title     
    h_fig = figure('Name', imageName);

    % Draw a plot showing the avg distance between nuclei per frame
    ymax = max (averageDist) + 1;
    plot (xAxis, averageDist); 
    title ('Average Distance between Cells (triangulated)');
    xlabel ('Frames');
    ylabel ('Avg distance (um)');
    if length (xAxis) > 1
       axis ([xAxis(1) xAxis(end) 0 ymax]);
    else
       axis ([xAxis(1) xAxis(1)+1 0 ymax]);
    end
end

% % Calculate the maximum distance value
% maxDist = max (max (distanceHist));
% 
% % Generate the bin vector by using the GUI field binsize
% bin = [];
% for binCount = 1 : binSize
%    binCol = round (binCount * maxDist / binSize);
%    bin = [bin binCol];
% end
% 
% % Generate a 3D distance histogram using 3-d bars which is chopped up in bins
% h_fig = figure ('Name', imageName);
% bar3 (bin, distanceHist, 0.5, 'detached');
% title ('3D Distance Histogram (binned)');