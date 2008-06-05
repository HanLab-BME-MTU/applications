function [pyxFiltered,vyxFiltered,vectorCoverageMask,runInfo]=removeVectorOutliers(runInfo,pyxVecKnown,vyxVecKnown)
% removeVectorOutliers removes outliers from a vector field

% runInfo is a structure containing at least these fields:
%
% thresh: search radius for distance matrix between all the vector origin
% pairs in raw data set. it should be big enough so that the average
% vector has several neighbors within thresh. thresh will be optimized
% based on vector density, so that on average there will be 7 neighbors
% for a given vector.
%
% nNeighborLowerLimit: vectors with a smaller number of neighbors will be
% considered outliers and removed from the filtered set. 3 is recommended.
%
% dirDev: mean directional deviation allowed (in degrees) for a vector and
% its neighbors. we find the mean(cos(theta)) for the vector and its
% neighbors and compare it to the corresponding value calculated from
% dirDev. vectors which deviate more than this value are considered
% outliers and removed from the filtered set.

global DEBUG__

% check input and assign parameters
if nargin<3
    error('Not enough input parameters')
end
if ~isstruct(runInfo)
    runInfo=struct;
end
if isfield(runInfo,'thresh')
    thresh=runInfo.thresh;
else
    thresh=30;
end
if isfield(runInfo,'nNeighborLowerLimit')
    nNeighborLowerLimit=runInfo.nNeighborLowerLimit;
else
    nNeighborLowerLimit=3;
end
if isfield(runInfo,'dirDev')
    dirDev=runInfo.dirDev;
else
    dirDev=10;
end


% get the distance between every vector pair in pyxVecKnown
D=createDistanceMatrix(pyxVecKnown,pyxVecKnown);

% inRange is 1 where vector (in y-direction) fits criteria
inRange=(D>0 & D<=thresh);
nNeighbors=sum(inRange,2);

% adapt thresh value until each vector has an average of 7 neighbors within
% thresh
meanNumNeighbors=mean(nNeighbors);
while meanNumNeighbors>7
    thresh=thresh-1;
    inRange=D>0 & D<=thresh;
    nNeighbors=sum(inRange,2);
    meanNumNeighbors=mean(nNeighbors);
end

% remove vectors without enough neighbors
inRange(nNeighbors<nNeighborLowerLimit,:)=0;
inRange=swapMaskValues(inRange,0,nan); % use NaN's instead of 0's

% using dot(A,B)=|A||B|cos(theta), we can find how coherent the
% vectors are.  coherency is measured by cos(theta) for each pair
mag=sqrt(sum(vyxVecKnown.^2,2)); % vector magnitudes
AB=mag*mag'; % matrix containing |A||B|
My=vyxVecKnown(:,1)*vyxVecKnown(:,1)'; % product of y-components
Mx=vyxVecKnown(:,2)*vyxVecKnown(:,2)'; % product of x-components
cosTheta=(My+Mx)./AB; 

% convert dirDev to radians and take cosine
cutOff=cos(dirDev*(pi()/180));

% remove incoherent vectors
inRange(nanmean(cosTheta.*inRange,2)<cutOff,:)=nan;

% retain any vector that still has a neighbor in inRange
idxFilt=nansum(inRange,2)>=1;
pyxFiltered=pyxVecKnown(idxFilt,:);
vyxFiltered=vyxVecKnown(idxFilt,:);

runInfo.meanNN=round(nanmean(nanmin(D.*inRange,[],2)));

% make mask with 1 on all pixels where a vector tail is located, then
% dilate/erode
vectorCoverageMask=zeros(runInfo.imL,runInfo.imW);
vectorCoverageMask(xy2index(pyxFiltered(:,2),pyxFiltered(:,1),runInfo.imL,runInfo.imW))=1;
maxNN=ceil(nanmax(nanmin(D.*inRange,[],2)));


vectorCoverageMask = imdilate(vectorCoverageMask,strel('disk',3*maxNN));
vectorCoverageMask = imerode(vectorCoverageMask,strel('disk',maxNN));
vectorCoverageMask=logical(vectorCoverageMask);

if DEBUG__==1
    % show mask negative with outliers in red and filtered set in yellow
    h=get(0,'CurrentFigure');
    if ~isempty(h);
        figure(h+1);
    else
        figure(1);
    end
    [M]=swapMaskValues(vectorCoverageMask,[0 1],[1 0]);
    imshow(M);
    hold on

    quiver(pyxVecKnown(:,2),pyxVecKnown(:,1),vyxVecKnown(:,2),vyxVecKnown(:,1),0,'r');
    quiver(pyxFiltered(:,2),pyxFiltered(:,1),vyxFiltered(:,2),vyxFiltered(:,1),0,'y');
end



