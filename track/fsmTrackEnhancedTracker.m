function currentM=fsmTrackEnhancedTracker(currentM,I,J,strg,counter1,gridSize,d0,userPath,threshold,influence,imgSize)
% fsmTrackMain uses the interpolated vector field to refine the tracking
%
% SYNOPSIS   currentM=fsmTrackEnhancedTracker(currentM,img,img2,strg,counter1,gridSize,d0,userPath,threshold,influence,imgSize)
%
% INPUT      currentM   : M matrix returned by the tracker for the current frame pair
%            I          : matrix [y x I]n of speckle coordinates and intensities for frame 1
%            J          : matrix [y x I]n of speckle coordinates and intensities for frame 2
%            strg       : format string for the correct numbering of files
%            counter1   : number of the first frame
%            gridSize   : distance between two interpolation points on the grid [gy gx]
%                         if gridSize is set equal to zero, the field is interpolated onto
%                         the original vector positions
%            d0         : correlation length or interpolation
%            userPath   : work directory (where the interpolated vector field is saved)
%            threshold  : radius of the region searched by the tracker for matching speckles
%            influence  : radius of influence for the neural network tracker. This is the initial search
%                         radius; the neural network can more effectively reconstruct flow if it can
%                         search over larger areas and therefore take more particles into account. The final
%                         matches will be constrained to have a maximum distance = 'radius', but initially 
%                         a larger search radius ('influence') will be used.
%            imgSize    : image size [y x]
%
% OUTPUT     currentM   : modified M matrix
%
% DEPENDENCES   fsmTrackEnhancedTracker uses { framework ; vectorFieldDiv ; updateD0FromDiv;
%                                              vectorFieldInterp ; fsmTrackPropSpecklePos ;
%                                              fsmTrackTrackerP }
%               fsmTrackEnhancedTracker is used by { fsmTrackMain }
%
% Aaron Ponti, January 14th, 2003

if nargin~=11
    error('Wrong number of input arguments');
end

% Extract vector field from currentM
raw=currentM(find(currentM(:,1)~=0 & currentM(:,3)~=0),:);

% In the very unlikely case that ALL particles are unmatched
if isempty(raw)
    
    % Create zero vectors
    vectors=[0 0 0 0];
    
    % Save interpolated vector field to disk for later use with gap closer
    indxStr=sprintf(strg,counter1);
    eval(['save ',userPath,filesep,'vectors',filesep,'vectors',indxStr,'.mat vectors;']);

    % Return the same matches
    return
    
end

if gridSize==0
    grid=raw(:,1:2); % Interpolate onto vector positions
else
    grid=framework(imgSize,gridSize); % Interpolate onto a regular grid
end

% Calculate the divergence of the vector field
[div,d0]=vectorFieldDiv(raw,grid,d0,[]);

% Update d0 depending on divergence
d0=updateD0FromDiv(div,d0,1,size(raw,1),size(grid,1));

% Interpolate vector field
vectors=vectorFieldInterp(raw,grid,d0,[]);

% Save interpolated vector field to disk for later use with gap closer
indxStr=sprintf(strg,counter1);
eval(['save ',userPath,filesep,'vectors',filesep,'vectors',indxStr,'.mat vectors;']);

% Propagate speckle positions
spPos=sortrows(currentM(find(currentM(:,1)~=0),1:2)); % Sorted positions
pSpPos=fsmTrackPropSpecklePos(spPos,vectors,'FORWARD');

% Create IP (with propagated positions)
Itmp=sortrows(I,1:2);
IP=cat(2,pSpPos,Itmp(:,3)); % Add original intensities to the propagated positions

% Track again with the brownian motion tracker with neural network
tmp2=fsmTrackTrackerBMTNN(IP,J,threshold,influence);

% Create a copy of tmp2
copyTmp2=tmp2;

% Correct back the coordinates of the first frame (= currentM(:,1:2))
for i=1:size(spPos,1)
    indx=find(tmp2(:,1)==pSpPos(i,1) & tmp2(:,2)==pSpPos(i,2));
    if length(indx)~=1
        error('Propagated position not univocally found in tmp2');
    end
    copyTmp2(indx,1:2)=spPos(i,:);
end
currentM=copyTmp2;
