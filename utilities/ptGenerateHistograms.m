function ptGenerateHistograms (ptPostpro, MPM, savePath, radioButtons)
% ptGenerateHistograms generates histograms for every frame of the movie. 
%
% SYNOPSIS       ptGenerateHistograms (ptPostpro, MPM)
%
% INPUT          ptPostpro : a structure which contains the information
%                            from the GUI
%                MPM       : matrix containing the cell tracks
%                savePath  : directory where the hist files will be saved
%                radioButtons : the values of the radiobuttons
%                
% OUTPUT         None (histograms are saved on disk) 
%
% DEPENDENCIES   ptGenerateHistograms   uses {nothing}
%                                  
%                ptsGenerateHistograms  is used by { PolyTrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Aug 04          Initial release
% Andre Kerstens        Sep 04          Added savePath as input parameter
% Andre Kerstens        Sep 04          Rewrite for new strcutures

% Get the latest data from the handles
cellProps = handles.allCellProps;
clusterProps = handles.allClusterProps;
frameProps = handles.allFrameProps;
jobData = handles.jobData;
guiData = handles.guiData;

% First assign all the postpro fields to a meaningfull variable
startFrame = handles.jobData;
endFrame = size(MPM/2,2);
increment = ptPostpro.increment;
plotStartFrame = ptPostpro.plotfirstimg;
plotEndFrame = ptPostpro.plotlastimg;
jobPath = ptPostpro.jobpath;
imageName = ptPostpro.imagenamenotiff;
numberOfFrames = ceil((plotEndFrame - plotStartFrame) / increment) + 1;
cellProps = ptPostpro.cellProps;
clusterProps = ptPostpro.clusterProps;
binSize = ptPostpro.binsize;
frameInterval = round (ptPostpro.timeperframe / 60);    % In minutes
pixelLength = ptPostpro.mmpixel;

% Get radio button values
generateAllCellsHist = ptPostpro.allcellshist;
generateSingleCellsHist = ptPostpro.singlecellshist;
generateClusteredCellsHist = ptPostpro.clusteredcellshist;

% In case the velocity over multiple frames is needed the following
% variable should be used
multipleFrameVelocity = ptPostpro.multframevelocity;

% Initialize the displacement and x-axis matrix
%displacementAll = zeros (size (MPM,1), numberOfFrames);
avgDisplacement = zeros (1, numberOfFrames-multipleFrameVelocity);
avgSingleDisplacement = zeros (1, numberOfFrames-multipleFrameVelocity);
avgClusteredDisplacement = zeros (1, numberOfFrames-multipleFrameVelocity);

% Initialize MPM counter
MPMCount = ceil ((plotStartFrame - startFrame) / increment) + 1;

% Initialize counter
iCount = 0;

% Initialize storage
avgDisplacement = [];
avgSingleDisplacement = [];
avgClusteredDisplacement = [];
avgVelocity = [];
avgSingleVelocity = [];
avgClusteredVelocity = [];
maxVelocity = [];
maxSingleCellVelocity = [];
maxClusteredCellVelocity = [];
histFrameCount = 10;

% Go through every frame of the set. Start at the second frame
% because only from there we can start calculating a displacement
for frameCount = plotStartFrame+increment : increment : plotEndFrame
    
   % Increase MPM counter
   MPMCount = MPMCount + 1;
   
   % If the velocity is calculated over multiple frames we should have
   % enough frames to do this
   if MPMCount > multipleFrameVelocity
   
      % Store the frame number for display on the x-axis
      iCount = iCount + 1;
   
      % Get the coordinates from the previous (dependent on the value of
      % multipleFrameVelocity) and the current frame
      prevCoordinates = MPM (:, 2*(MPMCount-multipleFrameVelocity)-1 : 2*(MPMCount-multipleFrameVelocity));
      curCoordinates = MPM (:, 2*MPMCount-1 : 2*MPMCount);
      coordinates = [prevCoordinates curCoordinates];
   
      % If coordinates in the current frame have started a track or the
      % ones in the previous frame ended a track, we cannot
      % calculate a displacement so we should throw these out
      coordinates (find ((coordinates (:,1) == 0 & coordinates (:,2) == 0) | ...
                         (coordinates (:,3) == 0 & coordinates (:,4) == 0)),:) = [];
      
      % From the cluster and cell properties we can find out which cells
      % are single and which belong to a cluster. First initialize two
      % matrices for this
      singleCellIndex = zeros (size (coordinates,1), 1);
      clusteredCellIndex = zeros (size (coordinates,1), 1);
      
      % Initialize single and clustered cell coordinates counters
      singleCellCount = 0;
      clusteredCellCount = 0;
      
      % Find the single cell and clustered coordinates and store these
      for jCount = 1 : size (coordinates,1) 

         % Find the index of the cell in CellProps
         cellIndex = find (cellProps (:,1,MPMCount) == coordinates (jCount,3) & ...
                           cellProps (:,2,MPMCount) == coordinates (jCount,4));
         
		 % In case the cell wasn't found in cellProps (maybe it was on the background?)
         % we can skip over this part
         if ~isempty (cellIndex)
             
            % Find this cluster in clusterProps
            clusterIndex = find (clusterProps (:,1,MPMCount) == cellProps (cellIndex,3,MPMCount));
         
            % Test whether this cell is in a cluster or not
            if clusterProps (clusterIndex,2,MPMCount) == 1  % Single cell found
            
               % Increase the counter for the single cell matrix
               singleCellCount = singleCellCount + 1; 
            
               % Store the coordinates
               singleCellIndex (singleCellCount,:) = jCount;
            else    % The cell was in a cluster
            
               % Increase the counter for the clustered cell matrix
               clusteredCellCount = clusteredCellCount + 1; 
            
               % Store the coordinates
               clusteredCellIndex (clusteredCellCount,:) = jCount;
            end
         end  % if ~isempty
      end  % for jCount = 1 : size (coordinates,1)
      
      % Remove all the zero entries in both matrices
      singleCellIndex (find (singleCellIndex(:,1) == 0), :) = [];
      clusteredCellIndex (find (clusteredCellIndex(:,1) == 0), :) = [];
      
      % Initialize displacement matrices
      displacement = zeros (size (coordinates,1),1);
      velocity = zeros (size (coordinates,1),1);
      singleCellsDisplacement = zeros (size (singleCellIndex,1),1);
      singleCellsVelocity = zeros (size (singleCellIndex,1),1);
      clusteredCellsDisplacement = zeros (size (clusteredCellIndex,1),1);
      clusteredCellsVelocity = zeros (size (clusteredCellIndex,1),1);  
      
      % Calculate the displacement for all cells
      displacement = sqrt ((coordinates(:,1) - coordinates(:,3)).^2 + ...
                           (coordinates(:,2) - coordinates(:,4)).^2);
                      
      % Calculate true velocity in um/min for all cells
      velocity = (displacement * pixelLength) / (multipleFrameVelocity * frameInterval);
             
      % Calculate true velocity in um/min for all cells
      singleCellVelocity = velocity(singleCellIndex,:);
            
      % Calculate true velocity in um/min for all cells
      clusteredCellVelocity = velocity(clusteredCellIndex);
            
      % Go to the save dir
      cd (savePath);
      
      % Calculate histogram for velocity all cells
      if generateAllCellsHist
         histVelAllCells = hist (velocity, binSize);
         filename = ['histVelAllCells_' num2str(iCount) '.mat'];
         save (filename, 'histVelAllCells');
      end
            
      % Save the velocity hist for single cells
      if generateSingleCellsHist
         histVelSingleCells = hist (singleCellVelocity, binSize);
         filename = ['histVelSingleCells_' num2str(iCount) '.mat'];
         save (filename, 'histVelSingleCells');
      end
      
      % Save the velocity hist for clustered cells
      if generateClusteredCellsHist
         histVelClusteredCells = hist (clusteredCellVelocity, binSize);
         filename = ['histVelClusteredCells_' num2str(iCount) '.mat'];
         save (filename, 'histVelClusteredCells');
      end
            
      % Store max velocity of the frame
      if ~isempty(velocity)
         maxVelocity(iCount) = max (velocity);
      else
         maxVelocity(iCount) = 0;
      end
      
      % Store max single cell velocity of the frame
      if ~isempty(singleCellVelocity)
         maxSingleCellVelocity(iCount) = max (singleCellVelocity);
      else
         maxSingleCellVelocity(iCount) = 0;
      end
      
      % Store max clustered cell velocity of the frame
      if ~isempty(clusteredCellVelocity)
         maxClusteredCellVelocity(iCount) = max (clusteredCellVelocity);
      else
         maxClusteredCellVelocity(iCount) = 0;
      end
      
   end   % if MPMCount > multipleFrameVelocity       
end   % for frameCount

% Save max velocity vectors
cd (savePath);
if generateAllCellsHist
   save ('maxVelAllCells.mat', 'maxVelocity');
end

if generateSingleCellsHist
   save ('maxVelSingleCells.mat', 'maxSingleCellVelocity');
end

if generateClusteredCellsHist
   save ('maxVelClusteredCells.mat', 'maxClusteredCellVelocity');
end
