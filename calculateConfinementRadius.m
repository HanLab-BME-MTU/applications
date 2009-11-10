function [experiment] = calculateConfinementRadius(experiment)
%calculateConfinementRadius outputs confinement radii for all tracks
%specified in input
%
%SYNOPSIS confinementRadius = calculateConfinementRadius(tracks,probDim);
%
%INPUT      experiment=   data structure pointing to all the
%                       image/detection/tracking data for a given condition;
%                       for e.g. 10 cell movies, the structure should have 10
%                       entries; the data structure can be created with the
%                       function loadConditionData, and it needs to contain
%                       at least the field
%                       .source, which is the (path) location of the
%                       lifetime information folder
%
%OUTPUT     confinementRadius:  Calculated the mean of the eigen vectors of
%                               the variance-covariance matrix of this
%                               track's positions
%           trackCenter:        calculate the track's center
%
%USES       none
%
%created by DA Nunez April 2, 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for endocytosis, always 2D
probDim = 2;
if nargin < 2 || isempty(rest)
   rest = [1 1 4 1 300]; 
end

%for each experiment
for iexp = 1:length(experiment)
    
    waitHandle = waitbar(iexp/length(experiment),['processing movie ' num2str(iexp) ' out of ' num2str(length(experiment))]);

    
    %Load Lifetime Information
    cd([experiment(iexp).source filesep 'LifetimeInfo'])
    load('lftInfo')
    % status matrix
    statMat = full(lftInfo.Mat_status);
    % lifetime matrix
    lftMat = full(lftInfo.Mat_lifetime);
    % x-coordinate matrix
    matX = full(lftInfo.Mat_xcoord);
    matX(matX == 0) =nan;
    matX(statMat == 4) = nan;
    matX(statMat == 5) = nan;
    % y-coordinate matrix
    matY = full(lftInfo.Mat_ycoord);
    matY(matY == 0) =nan;
    matY(statMat == 4) = nan;
    matY(statMat == 5) = nan;
    % disapp status matrix
    daMat = (lftInfo.Mat_disapp);
    % framerate
    framerate = experiment(iexp).framerate;

    %find all pits in movie that meet requirements specified by restriction
    %vector
    [findTrack,findFrame] = find((statMat==rest(1,1)) & (daMat==rest(1,2)) &...
        (lftMat>rest(1,3)) & (lftMat>round(rest(1,4)/framerate)) & (lftMat<round(rest(1,5)/framerate)));
    
    %estimate the confinement radius
    for iTrack = 1:size(findTrack,1)
        
        %get track coordinates
        xCoord = (matX(findTrack(iTrack),:))';
        yCoord = (matY(findTrack(iTrack),:))';
        xyzCoord = [xCoord yCoord];
        
        %find the eignevalues and eigenvectors of the variance-covariance
        %matrix of this track's positions
        eigenVal = eig(nancov(xyzCoord(:,1:probDim)));
        
        %calculate the track's confinement radius
        confRadius(iTrack,1) = sqrt( mean(eigenVal) * (probDim + 2) );
         if ~isreal(confRadius(iTrack,1))
             keyboard
         end
        
        %calculate the track's center
        trackCenter(iTrack,:) = nanmean(xyzCoord(:,1:probDim));
        
    end %for each track
    %store track confinement radii for each movie
    experiment(iexp).trackCenter = trackCenter;
    experiment(iexp).confRadius = confRadius;
    
    close(waitHandle)
end %of for each experiment

end %of function