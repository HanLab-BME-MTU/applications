%% Performs the actual motion clustering based on the estimated motion fields
% Output: filtered motion fields, clustering
% 


% function [] = clusterCellsCollectiveMotion(dirname,prefix,patchSize,nBilateralIter,t0,tn, padSize,colorRange,all)
function [] = clusterCellsCollectiveMotion(dirname,prefix,patchSize,nBilateralIter,t0,tn, padSize,colorRange, all)

% dirname = '/home/rack-wolf2/home/assafzar/WoundHealing/SN29/dic/11/auto';
% prefix = 'HGF_SN29_L11_Sum';

motionDir = strcat(dirname,'/MF/');
clustersDir = sprintf('%s/clusters/',motionDir);
clusterVisDir = sprintf('%s/vis/',clustersDir);

if (exist(sprintf('%s',clustersDir),'dir') == 0)
    eval(sprintf('mkdir %s',clustersDir));
end

if (exist(sprintf('%s',clusterVisDir),'dir') == 0)
    eval(sprintf('mkdir %s',clusterVisDir));
end

for i = t0 : tn-2
    clusterCells(dirname, sprintf('%s%s',prefix,pad(i,padSize)), patchSize, nBilateralIter, colorRange, all);
end



% load('/a/home/cc/cs/assafzar/clusteringMF/HGF_SN29_L11_Sum052_vf.mat');
% patchSize = 20;
% nIterations = 7;
% [filterMfsDxs,filterMfsDys] = bilateralFilteringDxy(dxs,dys,patchSize,nIterations);
% colorRange = [-10 10];
% figure; imagesc(filterMfsDxs); caxis(colorRange); title('Dx 7');
% figure; imagesc(filterMfsDys); caxis(colorRange); title('Dy 7');
end

% function [filterMfsDxs,filterMfsDys] = bilateralFilteringDxy(dxs,dys,patchSize,nIterations)
% [filterMfsDxs] = bilateralFiltering(dxs,patchSize,nIterations);
% [filterMfsDys] = bilateralFiltering(dys,patchSize,nIterations);
% 
% figure; imagesc(filterMfsDxs); title('Dx');
% figure; imagesc(filterMfsDys); title('Dy');
% 
% end

function [filterMfs] = bilateralFiltering(mfs,patchSize,nIterations)
maxVelocityOrig = max(abs(mfs(:)));
mfs =  mfs + maxVelocityOrig;
maxVelocityNew = max(mfs(:));
mfs =  mfs./maxVelocityNew;
lowResMfs = imresize(mfs,1.0/patchSize);
lowResMfs = max(lowResMfs,0);
lowResMfs = min(lowResMfs,1);

% bilateral filtering
filterMfs = lowResMfs;
for i = 1 : nIterations
    filterMfs = bfilter2(filterMfs,2);
end

filterMfs = filterMfs .* maxVelocityNew;
filterMfs = filterMfs - maxVelocityOrig;

% In case we want to go back to orinigal resolution
% newMfs = imresize(mfs,size(mfs));

end


%% 

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

%% Union of clusters in the same direction
% To unite to clusters the folowwing two condisions should be met:
% (1) less than 30/45/60 degrees between clusters
% % % (2) cells of one cluster "point" to the other cluster
% (2) close clusters
function [unitedClusters,unitedClustersMask] = unionClusters(clusters,dxs,dys)
clusterAngleTh = 30;
% clusterDirectionAngleTh = 30;
distanceTH = 50; % max distance between clusters to be united

unitedClusters.nclusters = 0;
unitedClusters.allClusters = {};
unitedClustersMask = zeros(size(dxs));
for cl1 = 1 : clusters.nclusters
    cluster1 = clusters.allClusters{cl1};
    cl1in = false;
    for cl2 = cl1+1 : clusters.nclusters        
        cluster2 = clusters.allClusters{cl2};
        
        dxsMask1 = dxs(cluster1);
        dysMask1 = dys(cluster1);
        dxsMask2 = dxs(cluster2);
        dysMask2 = dys(cluster2);
        
        v1 = [mean(dxsMask1(:)) mean(dysMask1(:))];
        v2 = [mean(dxsMask2(:)) mean(dysMask2(:))];
        
        % Measure angle between clusters and unify if < threshold
        
        % Normalize each vector bn = b/sqrt(sum(b.^2))
        vn1 = v1./sqrt(sum(v1.^2));
        vn2 = v2./sqrt(sum(v2.^2));
                                
        % Calc inner product: c = sum(an .* bn)
        innerProduct = sum(vn1 .* vn2);               
        % calc angle: acos(c) * 180/pi
        angle = acos(innerProduct) * 180/pi;
        
        % Unite!
        if angle < clusterAngleTh
            [clusterMask closecl] = closeClusters(cluster1,cluster2,distanceTH);            
            if closecl
                unitedClustersMask(clusterMask) = true;
                unitedClusters.nclusters = unitedClusters.nclusters + 1;
                unitedClusters.allClusters{unitedClusters.nclusters} = clusterMask;
                cl1in = true;
            end            
        end
    end
    if ~cl1in
        cluster1Erode = erode(cluster1,2); % erode so clusters will be seperated
        unitedClustersMask(cluster1Erode) = true; 
        unitedClusters.nclusters = unitedClusters.nclusters + 1;
        unitedClusters.allClusters{unitedClusters.nclusters} = cluster1Erode;
    end
end
end

%% Check whether the clusters are close
% closecl - true - close clusters
% clusterMask - merged clusters
function [clusterMask, closecl] = closeClusters(cluster1,cluster2,TH)
closecl = false;
clusterMask = false(size(cluster1));
DIST1 = bwdist(cluster1);
distC1C2 = DIST1;
distC1C2(~cluster2) = max(DIST1(:));
minDist = min(distC1C2(:)) + 1;
if minDist < TH    
    DIST2 = bwdist(cluster2);
    distC2C1 = DIST2;
    distC1C2(~cluster1) = max(DIST2(:));
    clusterMask = cluster1 | cluster2 | ((distC2C1 <= minDist) & (distC1C2 <= minDist)); 
    closecl = true;    
end
end

%%
function [] = clusterCells(dirname, prefix, patchSize, nIterations, colorRange, all)
% clear all;
close all;

imgDirname = sprintf('%s/images/',dirname);
dirname = sprintf('%s/MF/',dirname);

% dirname = '/home/rack-wolf2/home/assafzar/WoundHealing/SN29/dic/11/auto/MF/';
% patchSize = 20;
% nIterations = 7;
% colorRange = [-10 10];

clustersdir = sprintf('%s/clusters/',dirname);
visdir = sprintf('%s/clusters/vis/',dirname);

load(sprintf('%s/%s_vf.mat',dirname,prefix));

fbilateralFnameVisDx = sprintf('%s/%s_bilateralDX_%d.jpg',visdir,prefix,patchSize);
fbilateralFnameVisDy = sprintf('%s/%s_bilateralDY_%d.jpg',visdir,prefix,patchSize);

if ~all && exist(fbilateralFnameVisDx,'file')
    return;
end

[filterMfsDxs] = bilateralFiltering(dxs,patchSize,nIterations);
[filterMfsDys] = bilateralFiltering(dys,patchSize,nIterations);

h = figure; imagesc(filterMfsDxs); caxis(colorRange);
saveas(h,fbilateralFnameVisDx,'jpg');
close(h);

h = figure; imagesc(filterMfsDys); caxis(colorRange);
saveas(h,fbilateralFnameVisDy,'jpg');
close(h);

[outImgDx, outImgDy, unionFind] = RegionMerginSegmentationMF(filterMfsDxs,filterMfsDys);
outDxFname = sprintf('%s/%s_outSegDX_%d.jpg',visdir,prefix,patchSize);
h = figure; imagesc(outImgDx); caxis(colorRange);
saveas(h,outDxFname,'jpg');
close(h);
% imwrite(outImgDx,outDxFname,'jpg');
outDyFname = sprintf('%s/%s_outSegDY_%d.jpg',visdir,prefix,patchSize);
h = figure; imagesc(outImgDy); caxis(colorRange);
saveas(h,outDyFname,'jpg');
close(h);
% imwrite(outImgDy,outDyFname,'jpg');

motionMask = max(abs(outImgDx),abs(outImgDy)) > 2;
% fill holes and update union find
[parents] = parentsMap(unionFind,outImgDx);
[smallY,smallX] = size(parents);

% '/home/rack-wolf2/home/assafzar/WoundHealing/SN29/dic/11/auto/images/HGF_
% SN29_L11_Sum052.tif'
I = imread(sprintf('%s/%s.tif',imgDirname,prefix));


motionMask1 = imresize(motionMask,size(I),'nearest');
outImgDx1 = imresize(outImgDx,size(I),'nearest');
outImgDy1 = imresize(outImgDy,size(I),'nearest');
parents1 = imresize(parents,size(I),'nearest');

tic;
clusters.nclusters = 0;
clusters.allClusters = {};
ROI = false(size(parents1));
clustersMask = false(size(parents1));
for ind = 1 : smallY * smallX    
    cluster = parents1 == ind & motionMask1;    
    if sum(cluster(:)) > 8000   
        dxc = mean(outImgDx1(cluster));
        dyc = mean(outImgDy1(cluster));
        cluster = bwfill(cluster,'holes');        
        clusters.nclusters = clusters.nclusters + 1;        
        clusters.allClusters{clusters.nclusters} = cluster;
        ROI(cluster) = true;
        clustersMask(bwperim(cluster)) = true;
        % "fill holes" with correct values        
        outImgDx1(cluster) = dxc;
        outImgDy1(cluster) = dyc;        
    end
end
toc

ROI = imfill(ROI,'holes');

% [unitedClusters, unitedClustersMask] = unionClusters(clusters, outImgDx1, outImgDy1);

clustersInfoFname = sprintf('%s/%s_clustersInfo_%d.mat',clustersdir,prefix,patchSize);
save(clustersInfoFname,'filterMfsDxs','filterMfsDys', ...
    'outImgDx', 'outImgDy', 'unionFind', 'motionMask', 'parents', 'I', ...
    'outImgDx1', 'outImgDy1', 'clusters', 'ROI', 'clustersMask', 'parents1');%,'unitedClusters','unitedClustersMask');

Idisplay = I;
Idisplay(dilate(clustersMask,5)) = 255;

% '/a/home/cc/cs/assafzar/clusteringMF/HGF_SN29_L11_Sum052_clusters.jpg';
displayFname = sprintf('%s/%s_clusters_%d.jpg',visdir,prefix,patchSize);
jump = round(size(I,1)/30);
displayMotionFields(Idisplay,~ROI,outImgDy1,outImgDx1,jump,displayFname);

% original VF with clusters marked
% '/a/home/cc/cs/assafzar/clusteringMF/HGF_SN29_L11_Sum052_cinput.jpg';
displayInputFname = sprintf('%s/%s_cinput_%d.jpg',visdir,prefix,patchSize);
displayMotionFields(Idisplay,false(size(Idisplay)),dys,dxs,jump,displayInputFname);

% all cluster with vector fields
% '/a/home/cc/cs/assafzar/clusteringMF/HGF_SN29_L11_Sum052_allclusters.jpg';
IdisplayAll = I;
IdisplayAll(dilate(bwperim(ROI),5)) = 255;
displayAllFname = sprintf('%s/%s_allclusters_%d.jpg',visdir,prefix,patchSize);
displayMotionFields(IdisplayAll,~ROI,outImgDy1,outImgDx1,jump,displayAllFname);

% % all united cluster with vector fields
% IdisplayUnited = I;
% IdisplayUnited(dilate(bwperim(unitedClustersMask),5)) = 255;
% displayUnitedFname = sprintf('%s/%s_unitedclusters_%d.jpg',visdir,prefix,patchSize);
% displayMotionFields(IdisplayUnited,~ROI,outImgDy1,outImgDx1,jump,displayUnitedFname);
end