function [] = whCoordination(params,dirs)

time = 1 : params.nTime - params.frameJump;

fprintf('starting coordination\n');

for t = time   
    coordinationFname = [dirs.coordination pad(t,3) '_coordination.mat'];
    
    if exist(coordinationFname,'file') && ~params.always
        continue;
    end
    
    mfFname = [dirs.mfData pad(t,3) '_mf.mat']; % dxs, dys
        
    load(mfFname);
    
    sizeXY = size(dxs);
    
    [filterMfsDxs] = bilateralFiltering(dxs,params.patchSize,params.nBilateralIter);
    [filterMfsDys] = bilateralFiltering(dys,params.patchSize,params.nBilateralIter);

    %     dxs = imresize(dxs,1/params.patchSize);
    %     dys = imresize(dys,1/params.patchSize);
    
    [outImgDx, outImgDy, unionFind] = RegionMerginSegmentationMF(filterMfsDxs,filterMfsDys,params.regionMerginParams);
    [parents] = parentsMap(unionFind,outImgDx);
    parents(isnan(outImgDx)) = nan;
    [smallY,smallX] = size(parents);
    
    outImgDx1 = imresize(outImgDx,sizeXY,'nearest');
    outImgDy1 = imresize(outImgDy,sizeXY,'nearest');
    parents1 = imresize(parents,sizeXY,'nearest');    

    %% Filter small clusters
    clusters.nclusters = 0;
    clusters.allClusters = {};
    ROIclusters = false(size(parents1));
    clustersMask = false(size(parents1));
    for ind = 1 : smallY * smallX
        cluster = parents1 == ind;
        if sum(cluster(:)) > (params.minClusterSize / (params.pixelSize));
            dxc = mean(outImgDx1(cluster));
            dyc = mean(outImgDy1(cluster));
            cluster = bwfill(cluster,'holes');
            clusters.nclusters = clusters.nclusters + 1;
            clusters.allClusters{clusters.nclusters} = cluster;
            ROIclusters(cluster) = true;
            clustersMask(bwperim(cluster)) = true;
            % "fill holes" with correct values
            outImgDx1(cluster) = dxc;
            outImgDy1(cluster) = dyc;
        end
    end

    ROIclusters = imfill(ROIclusters,'holes');
    save(coordinationFname,'clusters','ROIclusters','clustersMask','outImgDx1','outImgDy1');
end

end


%% Bilateral filtering
function [filterMfs] = bilateralFiltering(mfs,patchSize,nIterations)
maxVelocityOrig = max(abs(mfs(:)));
mfs =  mfs + maxVelocityOrig;
maxVelocityNew = max(mfs(:));
mfs =  mfs./maxVelocityNew;
lowResMfs = imresize(mfs,1.0/patchSize);
nanMask = isnan(lowResMfs);
lowResMfs = max(lowResMfs,0);
lowResMfs = min(lowResMfs,1);

% bilateral filtering
filterMfs = lowResMfs;
for i = 1 : nIterations
    filterMfs = bfilter2(filterMfs,2);
end

filterMfs = filterMfs .* maxVelocityNew;
filterMfs = filterMfs - maxVelocityOrig;

filterMfs(nanMask) = nan;

% In case we want to go back to orinigal resolution
% newMfs = imresize(mfs,size(mfs));

end

%% For coordination
function [parents] = parentsMap(unionFind,I)
[sizeY,sizeX] = size(I);
parents = zeros(size(I));

for y = 1 : sizeY
    for x = 1 : sizeX
        parent = unionFind(y,x).parent;
        parents(y,x) = (parent.y-1) * sizeX + parent.x;
    end
end
end


%% get map to calculate strain rate
% xs = [1,  1,  0, -1,  -1,  -1,   0,   1];
% ys = [0,  1,  1,  1,   0,  -1,  -1,  -1];
%       0   45  90 -45   0    45  -90  -45
% >> atand(ys./xs)
% 
% ans =
% 
%      0    45    90   -45     0    45   -90   -45
function [strainDirectionsYs,strainDirectionsXs] = quantizeAngles(dys,dxs)
angles = atand(abs(dys./dxs)); 
quantAngles = 0:45:90;
nQuantAngles = length(quantAngles);

% Find closest angles
diff = nan(size(angles,1),size(angles,2),nQuantAngles);
for i = 1 : nQuantAngles
    diff(:,:,i) = abs(quantAngles(i) - angles);
end

[match,positions] = min(diff,[],3);

% % this is the closest quantified angle
% anglesQnatized = nan(size(angles));
% for i = 1 : nQuantAngles
%     anglesQnatized(positions == i) = quantAngles(i);    
% end

strainDirectionsXs = nan(size(angles));
strainDirectionsYs = nan(size(angles));

MASK = (positions == 1); % 0
strainDirectionsYs(MASK) = 0;
strainDirectionsXs(MASK & dxs >= 0) = 1;
strainDirectionsXs(MASK & dxs < 0) = -1;

MASK = (positions == 2); % 90
strainDirectionsXs(MASK & dxs >= 0) = 1;
strainDirectionsXs(MASK & dxs < 0) = -1;
strainDirectionsYs(MASK & dys >= 0) = 1;
strainDirectionsYs(MASK & dys < 0) = -1;

MASK = (positions == 3); % 90
strainDirectionsXs(MASK) = 0;
strainDirectionsYs(MASK & dys >= 0) = 1;
strainDirectionsYs(MASK & dys < 0) = -1;

% MASK = (positions == 1); % -90
% strainDirectionsXs(MASK) = 0;
% strainDirectionsYs(MASK) = -1;
% 
% MASK = (positions == 2); % -45
% strainDirectionsXs(MASK & dxs >= 0) = 1;
% strainDirectionsYs(MASK & dxs >= 0) = -1;
% strainDirectionsXs(MASK & dxs < 0) = -1;
% strainDirectionsYs(MASK & dxs < 0) = 1;
% 
% MASK = (positions == 3); % 0
% strainDirectionsYs(MASK) = 0;
% strainDirectionsXs(MASK & dxs >= 0) = 1;
% strainDirectionsXs(MASK & dxs < 0) = -1;
% 
% MASK = (positions == 4); % 45
% strainDirectionsXs(MASK & dxs >= 0) = 1;
% strainDirectionsYs(MASK & dxs >= 0) = 1;
% strainDirectionsXs(MASK & dxs < 0) = -1;
% strainDirectionsYs(MASK & dxs < 0) = -1;
% 
% MASK = (positions == 5); % 90
% strainDirectionsXs(MASK) = 0;
% strainDirectionsYs(MASK) = 1;
end

function [strainRate] = getStrainRate(directionX,directionY,speed)
[sizeY, sizeX] = size(speed);
[XGrid,YGrid] = meshgrid(1:size(speed,2),1:size(speed,1));

XGridBefore = XGrid -  directionX;
XGridAfter = XGrid +  directionX;
YGridBefore = YGrid -  directionY;
YGridAfter = YGrid +  directionY;
MASK = ...
    XGridBefore > 0 & XGridBefore < sizeX &...
    XGridAfter > 0 & XGridAfter < sizeX &...
    YGridBefore > 0 & YGridBefore < sizeY &...
    YGridAfter > 0 & YGridAfter < sizeY;

strainRate = nan(size(speed));
for y = 1 : sizeY
    for x = 1 : sizeX
        if MASK(y,x)
            strainRate(y,x) = speed(YGridAfter(y,x),XGridAfter(y,x)) - speed(YGridBefore(y,x),XGridBefore(y,x));
        end
    end
end

end



