function cenXYZ = centerMostPoint3D(maskIn,distX)
%CENTERMOSTPOINT3D finds the point in the input 3D mask which is furthest from its boundary
% 
% 
% coming soon... more documentation!
% 
% 
% 
% 
%Hunter Elliott 9/2011
%

%The centroid of this fraction of the highest distance-transform values
%will be used as the centermost point, as this is 
alpha = .01;
nPointsMin = 5;%Minimum number of points to average.

showPlots = true;

if nargin < 2 || isempty(distX)
    %Get distance transform, including distance to image boundary
    distX = bwdist(~padarray(maskIn,[1 1 1],0));
    distX = distX(2:end-1,2:end-1,2:end-1);
end


nPointsAvg = max(round(nnz(maskIn) * alpha),nPointsMin);

distThresh = sort(distX(maskIn(:)));

distThresh = distThresh(end-nPointsAvg);

centerCC = bwconncomp(distX >= distThresh);
areaSizes = cellfun(@numel,centerCC.PixelIdxList);
[~,iBiggest] = max(areaSizes);
[cenXYZ(:,2),cenXYZ(:,1),cenXYZ(:,3)] = ind2sub(size(maskIn),centerCC.PixelIdxList{iBiggest});

cenXYZ = mean(cenXYZ,1);

if showPlots
    figure;
    show3DMask(maskIn);
    hold on;
    spy3d(distX>=distThresh)
    plot3(cenXYZ(1),cenXYZ(2),cenXYZ(3),'or','MarkerSize',15)
end
    













