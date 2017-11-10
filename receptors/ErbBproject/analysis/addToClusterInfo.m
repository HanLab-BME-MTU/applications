function [clusterInfo] = addToClusterInfo(clusterInfo,points,varargin)
% addToClusterInfo add analysis to the results of meanshift clustering

ip = inputParser;
ip.CaseSensitive=false;
ip.StructExpand=true;

ip.addRequired('clusterInfo',@isstruct);
ip.addRequired('points',@(x) (isnumeric(x) & size(x,2) >= 2));

ip.addOptional('DoComposition',true,@islogical);
ip.addOptional('pixelSize',62.81,@isnumeric);

ip.parse(clusterInfo, points,varargin{:});

DoComposition = ip.Results.DoComposition;
pixelSize = ip.Results.pixelSize;

if DoComposition
    n = max(points(:,3));
end


for i=1:numel(clusterInfo)

    pnts = points(clusterInfo(i).ptIdData,:);

    if numel(unique(pnts(:,1))) > 2
    %finds convex hull and area
    [hull,area]= convhull(pnts(:,1:2));
    else
        hull =[];
        area = NaN;
    end
    
    % n is the number of pnts that are merged. Here we store the
    % number of points in this merged cluster from each element of the
    % shift array
    
    if DoComposition
        composition = hist(pnts(:,3),1:n);
        clusterInfo(i).composition = composition;
    end

    
    %Add values to clusterInfo
    clusterInfo(i).pnts = pnts;
    clusterInfo(i).area = area*pixelSize^2;
    clusterInfo(i).hull = hull;
    
end



end