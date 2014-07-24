function cenXYZ = centerMostPoint3D(maskIn,distX,varargin)
%CENTERMOSTPOINT3D finds the point in the input 3D mask which is furthest from its boundary
% 
% 
% coming soon... more documentation!
%
% Input:
%
%   maskIn - Mask with object to find centermost point int. Should contain
%            only one object (connected component).
%
%   distX - Distance transfrom. Optional.
% 
%   'PointType' - Specifies method to use for calculating centermost point.
%                 Either '3D' or '2.5D'. 
%                 3D uses full 3D mask and distance transform. 2D uses 2D
%                 projection of mask in Z direction and 2D distance
%                 transform to determine XY, and then uses max of distance
%                 along Z at this location to determine Z.
% 
%Hunter Elliott 9/2011
%

ip = inputParser;
ip.FunctionName = mfilename;
ip.addRequired('maskIn',@(x)(islogical(x) && ndims(x) == 3 && nnz(x) > 0));
ip.addOptional('distX',[],@(x)(ndims(x) == 3));
ip.addParamValue('FracUse',.01,@(x)(isscalar(x) && x >= 0 && x <= 1));%Fraction of centermost voxels to average for centermost point calc
ip.addParamValue('MinPointsUse',5,@(x)(isscalar(x) && isposint(x)));%Minimum number of points to use. Overrides the FracUse parameters
ip.addParamValue('PointType','3D',@(x)(ischar(x) && any(strcmpi(x,{'3D','2.5D'})))); %Type of centermost point calculation. 
%So input parser doesn't complain... stupid piece of shit
if ~isempty(varargin{:})
    ip.parse(maskIn,distX,varargin{:});
else
    ip.parse(maskIn,distX);
end
p = ip.Results;

showPlots = false;%For dev/test/debug

if nargin < 2 || isempty(distX)
    %Get distance transform, including distance to image boundary
    distX = bwdist(~padarray(maskIn,[1 1 1],0));
    distX = distX(2:end-1,2:end-1,2:end-1);
end


nPointsAvg = max(round(nnz(maskIn) * p.FracUse),p.MinPointsUse);

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
    












