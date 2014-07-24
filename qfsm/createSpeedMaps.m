function speedMap=createSpeedMaps(flow,n,corrLength,sampling,pixelSize,imgSize,gridSize,varargin)
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
% INPUT         flow      : the interpolated track matrix
%               n         : number of frames for temporal averaging (must be odd).
%                             If n==0, the whole stack is averaged into one speed map.
%               corrLength  : correlation length of the interpolator
%               sampling  : movie sampling size (s)
%               pixelSize : pixel size in the image domain (nm)
%               imgSize   : pixel size in the image domain (pixels)
%               gridSize
%               mask      : [ 0 | 1 ] overlays vector field to speed map
%
% OUTPUT        speedMap  : a cell array of matrices of size imgSize
%
%
% Aaron Ponti, January 16th, 2004
% Sebastien Besson, June 2011 (last modified March 2012)
% Adapted from fsmSpeedMaps


% Check input
ip = inputParser;
ip.addRequired('flow',@iscell);
ip.addRequired('n',@(x) isscalar(x) & mod(x,2)~=0);
ip.addRequired('corrLength',@isscalar);
ip.addRequired('sampling',@isscalar);
ip.addRequired('pixelSize',@isscalar);
ip.addRequired('imgSize',@isvector);
ip.addRequired('gridSize',@isscalar);
ip.addOptional('mask',true([imgSize numel(flow)+1]),@islogical);
ip.parse(flow,n,corrLength,sampling,pixelSize,imgSize,gridSize,varargin{:});
mask=ip.Results.mask;

% Get frame for interpolation
validFlowFrames = find(~cellfun(@isempty,flow));
assert(isequal(unique(diff(validFlowFrames)),1),...
    'Flow must be measured on consecutive frames to create speed maps');
startFrames = validFlowFrames(1:end-n+1);
midFrames = startFrames+fix(n/2);
preFrames=1:midFrames(1)-1;
postFrames = midFrames(end)+1:numel(flow)+1;

% Initialize  output
speedMap=cell(1,numel(flow)+1);


G=framework(imgSize,[gridSize gridSize]);
Md = zeros(size(G,1),4,numel(flow)+1);
for i=validFlowFrames
    Md(:,:,i) = vectorFieldAdaptInterp(flow{i},G,corrLength,[],'strain');
end
    
% Calculate average vector fields
for i=startFrames
    
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
speedMap(preFrames)=speedMap(midFrames(1));
speedMap(postFrames)=speedMap(midFrames(end));