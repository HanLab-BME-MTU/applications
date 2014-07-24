function idxPointsEdge = dissolveClustersCloseToEdge(obj,edgeDepth)

if edgeDepth > 0

% Determine slice size
nRandPoints = 1000;
z = repmat(mean(obj.data.points(:,3)),nRandPoints,1);
xMin = min(obj.data.points(:,1));
yMin = min(obj.data.points(:,2));
xRange = max(obj.data.points(:,1))-xMin;
yRange = max(obj.data.points(:,2))-yMin;
x = rand(nRandPoints,1)*xRange+xMin;
y = rand(nRandPoints,1)*yRange+yMin;
[~,dist] = KDTreeClosestPoint(obj.data.points,[x,y,z]);
sliceSize = mean(dist)^(1/3)*10;

% Density based noise filtering
filterSize = [sliceSize,sliceSize,2*(max(obj.data.points(:,3))-min(obj.data.points(:,3)))];
idxNeighborPoints = KDTreeRangeQuery(obj.data.points,obj.data.points,repmat(filterSize,obj.data.nPoints,1));
nNeighborPoints = cellfun(@numel,idxNeighborPoints);
nNeighborPointsAverage = mean(nNeighborPoints);
nNeighborPointThreshold = nNeighborPointsAverage*0.1; % 10 percent value is hardcoded
filteredPoints = obj.data.points(nNeighborPoints>nNeighborPointThreshold,:);

minX = min(filteredPoints(:,1)); maxX = max(filteredPoints(:,1));
minY = min(filteredPoints(:,2)); maxY = max(filteredPoints(:,2));
minZ = min(filteredPoints(:,3)); maxZ = max(filteredPoints(:,3));
extentX = maxX-minX; extentY = maxY-minY; extentZ = maxZ-minZ;

nSlicesX = ceil(extentX/sliceSize);
nSlicesY = ceil(extentY/sliceSize);

idxSlicesX = (1:nSlicesX)';
minXSlicesX = (idxSlicesX-1)*sliceSize+minX;
maxXSlicesX = idxSlicesX*sliceSize+minX;

idxSlicesY = (1:nSlicesY)';
minYSlicesY = (idxSlicesY-1)*sliceSize+minY;
maxYSlicesY = idxSlicesY*sliceSize+minY;

centersSlicesX = [(minXSlicesX+maxXSlicesX)/2,repmat((minY+maxY)/2,nSlicesX,1),repmat((minZ+maxZ)/2,nSlicesX,1)];
centersSlicesY = [repmat((minX+maxX)/2,nSlicesY,1),(minYSlicesY+maxYSlicesY)/2,repmat((minZ+maxZ)/2,nSlicesY,1)];

boxSizeSlicesX = repmat([sliceSize,extentY,extentZ],nSlicesX,1);
boxSizeSlicesY = repmat([extentX,sliceSize,extentZ],nSlicesY,1);

idxPointsInSlicesX = KDTreeRangeQuery(filteredPoints,centersSlicesX,boxSizeSlicesX);
idxPointsInSlicesY = KDTreeRangeQuery(filteredPoints,centersSlicesY,boxSizeSlicesY);

nPointsSlicesX = cellfun(@numel,idxPointsInSlicesX);
nPointsSlicesY = cellfun(@numel,idxPointsInSlicesY);
validSlicesX = nPointsSlicesX>0;
validSlicesY = nPointsSlicesY>0;
nValidSlicesX = nnz(validSlicesX);
nValidSlicesY = nnz(validSlicesY);

minYSlicesX = cellfun(@(a) min(filteredPoints(a,2)),idxPointsInSlicesX(validSlicesX));
maxYSlicesX = cellfun(@(a) max(filteredPoints(a,2)),idxPointsInSlicesX(validSlicesX));
minXSlicesY = cellfun(@(a) min(filteredPoints(a,1)),idxPointsInSlicesY(validSlicesY));
maxXSlicesY = cellfun(@(a) max(filteredPoints(a,1)),idxPointsInSlicesY(validSlicesY));

centersMinXSlicesY = [minXSlicesY+edgeDepth/2,centersSlicesY(validSlicesY,2:3)];
centersMaxXSlicesY = [maxXSlicesY-edgeDepth/2,centersSlicesY(validSlicesY,2:3)];
centersMinYSlicesX = [centersSlicesX(validSlicesX,1),minYSlicesX+edgeDepth/2,centersSlicesX(validSlicesX,3)];
centersMaxYSlicesX = [centersSlicesX(validSlicesX,1),maxYSlicesX-edgeDepth/2,centersSlicesX(validSlicesX,3)];
centersEdge = [centersMinXSlicesY;centersMaxXSlicesY;centersMinYSlicesX;centersMaxYSlicesX];

boxSizeMinXSlicesY = [repmat(edgeDepth,nValidSlicesY,1),boxSizeSlicesY(validSlicesY,2:3)];
boxSizeMaxXSlicesY = [repmat(edgeDepth,nValidSlicesY,1),boxSizeSlicesY(validSlicesY,2:3)];
boxSizeMinYSlicesX = [boxSizeSlicesX(validSlicesX,1),repmat(edgeDepth,nValidSlicesX,1),boxSizeSlicesX(validSlicesX,3)];
boxSizeMaxYSlicesX = [boxSizeSlicesX(validSlicesX,1),repmat(edgeDepth,nValidSlicesX,1),boxSizeSlicesX(validSlicesX,3)];
boxSizeEdge = [boxSizeMinXSlicesY;boxSizeMaxXSlicesY;boxSizeMinYSlicesX;boxSizeMaxYSlicesX];

% slice = 10;
% plot(filteredPoints(:,1),filteredPoints(:,2), 'b.');
% hold on;
% plot(centersSlicesX(slice,1) + boxSizeSlicesX(slice,1)/2*[-1 -1 1 1 -1],centersSlicesX(slice,2) + boxSizeSlicesX(slice,2)/2*[1 -1 -1 1 1], 'g-');
% plot(filteredPoints(idxPointsInSlicesX{slice},1), filteredPoints(idxPointsInSlicesX{slice},2), 'r.');
% for i=1:nSlicesY
%     slice = i;
%     plot(centersMinXSlicesY(slice,1) + boxSizeMinXSlicesY(slice,1)/2*[-1 -1 1 1 -1],centersMinXSlicesY(slice,2) + boxSizeMinXSlicesY(slice,2)/2*[1 -1 -1 1 1], 'k-');
%     plot(centersMaxXSlicesY(slice,1) + boxSizeMaxXSlicesY(slice,1)/2*[-1 -1 1 1 -1],centersMaxXSlicesY(slice,2) + boxSizeMaxXSlicesY(slice,2)/2*[1 -1 -1 1 1], 'k-');
% end
% for i=1:nSlicesX
%     slice = i;
%     plot(centersMinYSlicesX(slice,1) + boxSizeMinYSlicesX(slice,1)/2*[-1 -1 1 1 -1],centersMinYSlicesX(slice,2) + boxSizeMinYSlicesX(slice,2)/2*[1 -1 -1 1 1], 'k-');
%     plot(centersMaxYSlicesX(slice,1) + boxSizeMaxYSlicesX(slice,1)/2*[-1 -1 1 1 -1],centersMaxYSlicesX(slice,2) + boxSizeMaxYSlicesX(slice,2)/2*[1 -1 -1 1 1], 'k-');
% end
% axis equal

idxPointsEdge = KDTreeRangeQuery(obj.data.points,centersEdge,boxSizeEdge);
idxPointsEdge = vertcat(idxPointsEdge{:});

% Append the low density points
% idxPointsEdge = [idxPointsEdge;find(nNeighborPoints<=nNeighborPointThreshold)];

idxPointsEdge = unique(idxPointsEdge);

% Dissolve clusters containing points with an index in idxPointsEdge
nEdgePointsInCluster = cellfun(@(a) nnz(ismember(a,idxPointsEdge)),obj.data.clusters);

cluster2Keep = obj.data.clusters(nEdgePointsInCluster==0);
cluster2Dissolve = obj.data.clusters(nEdgePointsInCluster>0);
cluster2Dissolve = num2cell(horzcat(cluster2Dissolve{:})');
obj.data.clusters = [cluster2Keep;cluster2Dissolve];

modelType2Keep = obj.data.modelType(nEdgePointsInCluster==0);
modelType2Dissolve = zeros(numel(cluster2Dissolve),1);
obj.data.modelType = [modelType2Keep;modelType2Dissolve];

modelLength2Keep = obj.data.modelLength(nEdgePointsInCluster==0);
modelLength2Dissolve = zeros(numel(cluster2Dissolve),1);
obj.data.modelLength = [modelLength2Keep;modelLength2Dissolve];

modelRes2Keep = obj.data.modelRes(nEdgePointsInCluster==0);
modelRes2Dissolve = cell(numel(cluster2Dissolve),1);
obj.data.modelRes = [modelRes2Keep;modelRes2Dissolve];

modelProj2Keep = obj.data.modelProj(nEdgePointsInCluster==0);
modelProj2Dissolve = cell(numel(cluster2Dissolve),1);
obj.data.modelProj = [modelProj2Keep;modelProj2Dissolve];

modelBezCP2Keep = obj.data.modelBezCP(nEdgePointsInCluster==0);
modelBezCP2Dissolve = cell(numel(cluster2Dissolve),1);
obj.data.modelBezCP = [modelBezCP2Keep;modelBezCP2Dissolve];

clusterColor2Keep = obj.data.clusterColor(nEdgePointsInCluster==0,:);
clusterColor2Dissolve = repmat(obj.data.nullClusterColor,numel(cluster2Dissolve),1);
obj.data.clusterColor = [clusterColor2Keep;clusterColor2Dissolve];

obj.updateClusterData();

disp('GeometricMatcher: Clusters close to an edge dissolved!');
end

end

