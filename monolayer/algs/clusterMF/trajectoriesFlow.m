% Simulate repeatedly trajectories and display flows
function [] = trajectoriesFlow(dirname,prefix,patchSize)

motionDir = strcat(dirname,'/MF/');
trajClustersDir = sprintf('%s/trajClusters/',motionDir);
trajFlowDir = sprintf('%s/trajFlow/',trajClustersDir);

if (exist(sprintf('%s',trajFlowDir),'dir') == 0)
    eval(sprintf('mkdir %s',trajFlowDir));
end

simulateFlow(dirname, prefix, 'Phase1', patchSize);
simulateFlow(dirname, prefix, 'Phase2', patchSize);
simulateFlow(dirname, prefix, 'Phase3', patchSize);
end

%%
function [] = simulateFlow(dirname, prefix, phaseStr, patchSize)

% clear all;
close all;

motionDir = sprintf('%s/MF/',dirname);
trajClustersDir = sprintf('%s/trajClusters/',motionDir);
trajFlowDir = sprintf('%s/trajFlow/',trajClustersDir);

clustersInfoFname = sprintf('%s/clusters/%s_clustersTrajInfo_%s.mat',trajClustersDir,prefix,phaseStr);

load(clustersInfoFname);
% 'dxs1','dys1','outImgDxTraj', 'outImgDyTraj', 'unionFindTraj', 'motionMaskTraj', 'parentsTraj', 'It0', 'I', ...
% 'outImgDx1Traj', 'outImgDy1Traj', 'clusters', 'ROI', 'clustersMask', 'parents1Traj',...
% 'ROI_Seg','ROI_Seg1'

SIMULATION_ITERS = 100;
[sizeY, sizeX] = size(I);
rangeY = round(patchSize/2) : patchSize : sizeY;
rangeX = round(patchSize/2) : patchSize : sizeX;
startYs = repmat(rangeY,length(rangeX),1)';
startXs = repmat(rangeX,length(rangeY),1);
endYs = startYs;
endXs = startXs;

% roiFname = sprintf('%s/ROI/%s000_ROI.mat',dirname,prefix);
% load(roiFname);
% ROI = binNoExcesses;
% [cogX, cogY] = cog(ROI);
% deadZone = false(size(ROI));
% deadZone(cogY,:) = true;
% deadZone = dilate(deadZone,150);

for iter = 1 : SIMULATION_ITERS
    for y = 1 : size(endYs,1)
        for x = 1 : size(endYs,2)
            if ROI_Seg(y,x)
                curY = endYs(y,x);
                curX = endXs(y,x);
%                 if ~deadZone(round(curY),round(curX))
                    endYs(y,x) = curY + dys1(round(curY),round(curX));
                    endXs(y,x) = curX + dxs1(round(curY),round(curX));
%                 end
            end
        end
    end
    endYs = min(max(endYs,1),sizeY);
    endXs = min(max(endXs,1),sizeX);
end


% %% clustering
% [outImgDxTraj, outImgDyTraj, unionFindTraj] = RegionMerginSegmentationMF(dxs,dys);

finalLocationMapFname = sprintf('%s/%s_FlowMap_%s.jpg',trajFlowDir,prefix,phaseStr);
endXs1 = imresize(endXs,size(I),'nearest');
h = figure; imagesc(endXs1); axis([0,1024,0,1024]); colorbar; 
saveas(h,finalLocationMapFname,'jpg');
close(h);

finalLocationMapYsFname = sprintf('%s/%s_FlowYs_%s.jpg',trajFlowDir,prefix,phaseStr);
endYs1 = imresize(endYs,size(I),'nearest');
h = figure; imagesc(endYs1); axis([0,1024,0,1024]); colorbar; 
saveas(h,finalLocationMapYsFname,'jpg');
close(h);

% Trajectories estimation visualization
dxs = endXs - startXs;
dys = endYs - startYs;
dxs1 = imresize(dxs,size(I),'nearest');
dys1 = imresize(dys,size(I),'nearest');
displayFname = sprintf('%s/%s_FlowVis_%s.jpg',trajFlowDir,prefix,phaseStr);
jump = round(size(I,1)/30);
displayMotionFields(It0,false(size(ROI)),dys1,dxs1,jump,displayFname);

%% cluster flow
params.P = 0.03;
params.Q = 0.001; %small Q --> less merging
[outX, outY, uf] = RegionMerginSegmentationMF(endXs,zeros(size(endXs)),params);

finalLocationMapClustersFname = sprintf('%s/%s_FlowMapClusters_%s.jpg',trajFlowDir,prefix,phaseStr);
outX1 = imresize(outX,size(I),'nearest');
h = figure; imagesc(outX1); axis([0,1024,0,1024]); colorbar; 
saveas(h,finalLocationMapClustersFname,'jpg');
close(h);

% clusters.nclusters = 0;
% clusters.allClusters = {};
% ROI = false(size(parents1Traj));
% clustersMask = false(size(parents1Traj));
% for ind = 1 : smallY * smallX    
%     cluster = parents1Traj == ind;    
%     if sum(cluster(:)) > 4000   
%         dxc = mean(outImgDx1Traj(cluster));
%         dyc = mean(outImgDy1Traj(cluster));
%         cluster = bwfill(cluster,'holes');        
%         clusters.nclusters = clusters.nclusters + 1;        
%         clusters.allClusters{clusters.nclusters} = cluster;
%         ROI(cluster) = true;
%         clustersMask(bwperim(cluster)) = true;
%         % "fill holes" with correct values        
%         outImgDx1Traj(cluster) = dxc;
%         outImgDy1Traj(cluster) = dyc;        
%     end
% end
% 
% IFlowClusters = It0;


end