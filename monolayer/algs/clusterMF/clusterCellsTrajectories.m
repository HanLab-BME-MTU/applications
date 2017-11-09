% Performs the clustering based on cells' trajectories
function [] = clusterCellsTrajectories(dirname,prefix,patchSize,t0, phase1, phase2, phase3, padSize, colorRange, all)
% dirname = '/home/rack-wolf2/home/assafzar/WoundHealing/SN29/dic/11/auto';
% prefix = 'HGF_SN29_L11_Sum';

motionDir = strcat(dirname,'/MF/');
clustersDir = sprintf('%s/clusters/',motionDir);
clusterVisDir = sprintf('%s/vis/',clustersDir);
trajClustersDir = sprintf('%s/trajClusters/',motionDir);
trajClustersDirVis = sprintf('%s/trajClusters/vis/',motionDir);
trajClustersDirClusters = sprintf('%s/trajClusters/clusters/',motionDir);

if (exist(sprintf('%s',clustersDir),'dir') == 0)
    eval(sprintf('mkdir %s',clustersDir));
end

if (exist(sprintf('%s',clusterVisDir),'dir') == 0)
    eval(sprintf('mkdir %s',clusterVisDir));
end

if (exist(sprintf('%s',trajClustersDir),'dir') == 0)
    eval(sprintf('mkdir %s',trajClustersDir));
    eval(sprintf('mkdir %s',trajClustersDirVis));
    eval(sprintf('mkdir %s',trajClustersDirClusters));
end

clusterTrajectories(dirname, prefix, 'Phase1', t0, phase1, patchSize, padSize, colorRange, all);
clusterTrajectories(dirname, prefix, 'Phase2', phase1+1, phase2, patchSize, padSize, colorRange, all);
clusterTrajectories(dirname, prefix, 'Phase3', phase2+1, phase3-2, patchSize, padSize, colorRange, all);
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
function [] = clusterTrajectories(dirname, prefix, phaseStr, tstart, tend, patchSize, padSize, colorRange, all)
% clear all;
close all;
all = false;

imgDirname = sprintf('%s/images/',dirname);
motionDir = sprintf('%s/MF/',dirname);

bilateralFilterDir = sprintf('%s/clusters/',motionDir);
trajectoriesCustersDir = sprintf('%s/trajClusters/',motionDir);
trajClustersDirClusters = sprintf('%s/clusters/',trajectoriesCustersDir);
visdir = sprintf('%s/vis/',trajectoriesCustersDir);

% initiate
load(sprintf('%s/%s%s_clustersInfo_%d.mat',bilateralFilterDir,prefix,pad(tstart,padSize),patchSize));
It0 = I;
[sizeY, sizeX] = size(I);
rangeY = round(patchSize/2) : patchSize : sizeY;
rangeX = round(patchSize/2) : patchSize : sizeX;
startTrajYs = repmat(rangeY,length(rangeX),1)';
startTrajXs = repmat(rangeX,length(rangeY),1);
endTrajYs = startTrajYs;
endTrajXs = startTrajXs;

outDxFname = sprintf('%s/%s_outTrajSegDX_%s.jpg',visdir,prefix,phaseStr);
outDyFname = sprintf('%s/%s_outTrajSegDY_%s.jpg',visdir,prefix,phaseStr);
inDxFname = sprintf('%s/%s_inTrajSegDX_%s.jpg',visdir,prefix,phaseStr);
inDyFname = sprintf('%s/%s_inTrajSegDY_%s.jpg',visdir,prefix,phaseStr);

% ROI
% roiFname = sprintf('%s/ROI/HGF_SN29_L11_Sum%s_ROI.mat',dirname,pad(tstart,padSize)); % FOR TESTING
roiFname = sprintf('%s/ROI/%s%s_ROI.mat',dirname,prefix,pad(tstart,padSize));
if (exist(roiFname,'file') == 0)
    ROI_Seg1 = true(size(I));
else
    load(roiFname);
    ROI_Seg1 = ~binNoExcesses;
    clear binNoExcesses;
end
ROI_Seg = imresize(ROI_Seg1,1.0/15);
ROI_Seg = ROI_Seg(1:end-1,1:end-1);

clustersInfoFname = sprintf('%s/%s_clustersTrajInfo_%s.mat',trajClustersDirClusters,prefix,phaseStr);

if all || exist(clustersInfoFname,'file') == 0
    for t = tstart : tend
        bilateralFilterFname = sprintf('%s/%s%s_clustersInfo_%d.mat',bilateralFilterDir,prefix,pad(t,padSize),patchSize);
        load(bilateralFilterFname); % outImgDx1, outImgDy1
        for y = 1 : size(endTrajYs,1)
            for x = 1 : size(endTrajYs,2)
                if ROI_Seg(y,x)
                    curY = endTrajYs(y,x);
                    curX = endTrajXs(y,x);
                    endTrajYs(y,x) = curY + outImgDy1(round(curY),round(curX));
                    endTrajXs(y,x) = curX + outImgDx1(round(curY),round(curX));
                end
            end
        end
        endTrajYs = min(max(endTrajYs,1),sizeY);
        endTrajXs = min(max(endTrajXs,1),sizeX);
    end
else
    load(outDxFname);
    load(outDyFname);
end

dxs = endTrajXs - startTrajXs;
dys = endTrajYs - startTrajYs;
%% clustering
[outImgDxTraj, outImgDyTraj, unionFindTraj] = RegionMerginSegmentationMF(dxs,dys);

%%
motionMaskTraj = abs(outImgDxTraj) > 1*(tend-tstart) | abs(outImgDyTraj) > 1*(tend-tstart);
% fill holes and update union find
[parentsTraj] = parentsMap(unionFindTraj,outImgDxTraj);
[smallY,smallX] = size(parents);

motionMask1Traj = imresize(motionMaskTraj,size(I),'nearest');
outImgDx1Traj = imresize(outImgDxTraj,size(I),'nearest');
outImgDy1Traj = imresize(outImgDyTraj,size(I),'nearest');
parents1Traj = imresize(parentsTraj,size(I),'nearest');
dxs1 = imresize(dxs,size(I),'nearest');
dys1 = imresize(dys,size(I),'nearest');

%%
h = figure; imagesc(dxs); colorbar; caxis([-150,150].*1.24);%caxis(colorRange);
saveas(h,inDxFname,'jpg');
close(h);

h = figure; imagesc(dys); colorbar; caxis([-150,150]*1.24);%caxis(colorRange);
saveas(h,inDyFname,'jpg');
close(h);


h = figure; imagesc(outImgDxTraj); colorbar; caxis([-150,150]*1.24);%caxis(colorRange);
saveas(h,outDxFname,'jpg');
close(h);

h = figure; imagesc(outImgDyTraj); colorbar; caxis([-150,150]*1.24);%caxis(colorRange);
saveas(h,outDyFname,'jpg');
close(h);

%% Filter small clusters

clusters.nclusters = 0;
clusters.allClusters = {};
ROI = false(size(parents1Traj));
clustersMask = false(size(parents1Traj));
for ind = 1 : smallY * smallX    
    cluster = parents1Traj == ind & motionMask1Traj;    
    if sum(cluster(:)) > 4000   
        dxc = mean(outImgDx1Traj(cluster));
        dyc = mean(outImgDy1Traj(cluster));
        cluster = bwfill(cluster,'holes');        
        clusters.nclusters = clusters.nclusters + 1;        
        clusters.allClusters{clusters.nclusters} = cluster;
        ROI(cluster) = true;
        clustersMask(bwperim(cluster)) = true;
        % "fill holes" with correct values        
        outImgDx1Traj(cluster) = dxc;
        outImgDy1Traj(cluster) = dyc;        
    end
end

ROI = imfill(ROI,'holes');

% [unitedClusters, unitedClustersMask] = unionClusters(clusters, outImgDx1, outImgDy1);

save(clustersInfoFname,'dxs1','dys1','outImgDxTraj', 'outImgDyTraj', 'unionFindTraj', 'motionMaskTraj', 'parentsTraj', 'It0', 'I', ...
    'outImgDx1Traj', 'outImgDy1Traj', 'clusters', 'ROI', 'clustersMask', 'parents1Traj',...
    'ROI_Seg','ROI_Seg1');

Idisplay = It0;
Idisplay(dilate(clustersMask,5)) = 255;

% Trajectories estimation visualization
displayFname = sprintf('%s/%s_trajectories_%s.jpg',visdir,prefix,phaseStr);
jump = round(size(I,1)/30);
displayMotionFields(It0,false(size(ROI)),dys1,dxs1,jump,displayFname);

% clusters over original motion fields
displayFname = sprintf('%s/%s_initclustersTraj_%s.jpg',visdir,prefix,phaseStr);
jump = round(size(I,1)/30);
displayMotionFields(Idisplay,false(size(ROI)),dys1,dxs1,jump,displayFname);

% '/a/home/cc/cs/assafzar/clusteringMF/HGF_SN29_L11_Sum052_clusters.jpg';
displayFname = sprintf('%s/%s_clustersTraj_%s.jpg',visdir,prefix,phaseStr);
jump = round(size(I,1)/30);
displayMotionFields(Idisplay,~ROI,outImgDy1Traj,outImgDx1Traj,jump,displayFname);

% all cluster with vector fields
% '/a/home/cc/cs/assafzar/clusteringMF/HGF_SN29_L11_Sum052_allclusters.jpg';
IdisplayAll = I;
IdisplayAll(dilate(bwperim(ROI),5)) = 255;
displayAllFname = sprintf('%s/%s_allclustersTraj_%s.jpg',visdir,prefix,phaseStr);
displayMotionFields(IdisplayAll,~ROI,outImgDy1Traj,outImgDx1Traj,jump,displayAllFname);

% % all united cluster with vector fields
% IdisplayUnited = I;
% IdisplayUnited(dilate(bwperim(unitedClustersMask),5)) = 255;
% displayUnitedFname = sprintf('%s/%s_unitedclusters_%d.jpg',visdir,prefix,patchSize);
% displayMotionFields(IdisplayUnited,~ROI,outImgDy1,outImgDx1,jump,displayUnitedFname);

%% Circular statistics
XMAG = [];
YMAG = [];
for y = 1 : size(endTrajYs,1)
    for x = 1 : size(endTrajYs,2)
        if ROI_Seg(y,x) && motionMaskTraj(y,x)
            XMAG = [XMAG dxs(y,x)];
            YMAG = [YMAG dys(y,x)];            
        end
    end
end
RES = xy2rad([XMAG',YMAG']);
RADIANS = RES(:,1); % [angle,mag,xcor,ycor]
circHistBins = 40;
[h,circHist] = circ_plot(RADIANS, 'hist',[],circHistBins,true,false);%format, formats, varargin

circStatsIFname = sprintf('%s/%s_circStats_%s.jpg',visdir,prefix,phaseStr);
saveas(h,circStatsIFname,'jpg');
circStatsFname = sprintf('%s/%s_circStats_%s.mat',trajClustersDirClusters,prefix,phaseStr);
save(circStatsFname,'XMAG','YMAG','RADIANS','circHist');
end