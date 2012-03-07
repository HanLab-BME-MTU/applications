function speedMap=createSpeedMaps(Md,n,sampling,pixelSize,imgSize,mask)
% fsmSpeedMaps creates speed maps from the flow maps
%
% fsmSpeedMaps goes through the whole M (or Md) stack of s matrices (each 
%  matrix corresponds to the    matches returned by the tracker for frames
%  2 consecutive frames) and creates t<=s speed maps  each of which is the
%  average over n frames [j-n/2:j+n/2] around frame j, j in 
%  [fix(n/2)+1:s-fix(n/2)]
%   
% SYNOPSIS      outputdir=fsmSpeedMaps(Md,n,sampling,pixelSize,imgSize,mask)
%
% INPUT         Md        : the interpolated track matrix
%               n         : number of frames for temporal averaging (must be odd).
%                             If n==0, the whole stack is averaged into one speed map.
%               sampling  : movie sampling size (s)
%               pixelSize : pixel size in the image domain (nm)
%               imgSize   : pixel size in the image domain (pixels)
%               mask      : [ 0 | 1 ] overlays vector field to speed map
%
% OUTPUT        speedMap  : a cell array of matrices of size imgSize
%
%
% Aaron Ponti, January 16th, 2004
% Sebastien Besson, June 2011
% Adapted from fsmSpeedMaps

% Calculate average vector fields
nMaps=size(Md,3)-n+1;
speedMap=cell(1,size(Md,3)+1);
for i=1:nMaps
    
    % Extract mean velocities over n frames
    vectors=Md(:,3:4,i:i+n-1)-Md(:,1:2,i:i+n-1);
    meanVectors=mean(vectors,3);
    
    % Construct average vector matrix
    Mav=zeros(size(Md,1),size(Md,2));
    Mav(:,1:2)=Md(:,1:2,1);
    Mav(:,3:4)=Mav(:,1:2)+meanVectors;
    
    % Calculate velocities
    velocities=zeros(size(Mav,1),3);
    velocities(:,1:2)=Mav(:,1:2);
    velocities(:,3)=sqrt(meanVectors(:,1).^2+meanVectors(:,2).^2);
    
    % Reshape
    iMap=i+fix(n/2);
    speedMap{iMap}=reshape(velocities(:,3),length(unique(Md(:,2,1))),length(unique(Md(:,1,1))))';
    
    % Transform to nm/min
    speedMap{iMap}=speedMap{iMap}*(60/sampling)*pixelSize;
    
    % If needed, resize
    speedMap{iMap}=imresize(speedMap{iMap},imgSize,'bilinear');

    % Use union of masks
    speedMap{iMap} = speedMap{iMap} .* logical(sum(mask(:,:,i:i+n-1),3));
end
speedMap(1:fix(n/2))=speedMap(fix(n/2)+1);
speedMap(nMaps+fix(n/2)+1:end)=speedMap(nMaps+fix(n/2));