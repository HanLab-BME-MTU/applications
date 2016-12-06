function [outImgDx, outImgDy, unionFind] = RegionMerginSegmentationMF(dxs,dys,params)

if nargin < 3
    params.P = 0.03;% small P --> more merging
    params.Q = 0.01;% large Q --> more merging (more significant than P)
end

% global P;
% global Q;
global contLP;
global logTable;

[sizeY,sizeX] = size(dxs);
outImgDx = zeros(sizeY,sizeX);
outImgDy = zeros(sizeY,sizeX);

% P = 0.03; % small P --> more merging
% Q = 0.01; % large Q --> more merging (more significant than P)
contLP = log(2.0/params.P);
logTable = log((1:(sizeY*sizeX) )+ 1);

% Initiate all pairs
numGroups = 2*(sizeY-2) * (sizeX-2);
allPairs(numGroups).x1 = 0;
allPairs(numGroups).y1 = 0;
allPairs(numGroups).x2 = 0;
allPairs(numGroups).y2 = 0;
allPairs(numGroups).diffDx = 0;
allPairs(numGroups).diffDy = 0;
allPairs(numGroups).diff = 0;
i = 1;
for y = 2 : sizeY - 1
    for x = 2 : sizeX - 1  
        % right
        allPairs(i).x1 = x;
        allPairs(i).y1 = y;
        allPairs(i).x2 = x+1;
        allPairs(i).y2 = y;        
        allPairs(i).diffDx = abs(dxs(y,x) - dxs(y,x+1)) / max(abs(dxs(y,x)),abs(dxs(y,x+1)));
        allPairs(i).diffDy = abs(dys(y,x) - dys(y,x+1)) / max(abs(dys(y,x)),abs(dys(y,x+1)));
%         allPairs(i).diff = max(allPairs(i).diffDx,allPairs(i).diffDy);
        allPairs(i).diff = mean([allPairs(i).diffDx,allPairs(i).diffDy]);
        i = i + 1; % next group
        
        % down
        allPairs(i).x1 = x;
        allPairs(i).y1 = y;
        allPairs(i).x2 = x;
        allPairs(i).y2 = y+1;
        allPairs(i).diffDx = abs(dxs(y,x) - dxs(y+1,x)) / max(abs(dxs(y,x)),abs(dxs(y+1,x)));
        allPairs(i).diffDy = abs(dys(y,x) - dys(y+1,x)) / max(abs(dys(y,x)),abs(dys(y+1,x)));
%         allPairs(i).diff = max(allPairs(i).diffDx,allPairs(i).diffDy);
        allPairs(i).diff = mean([allPairs(i).diffDx,allPairs(i).diffDy]);
        i = i + 1; % next group
    end
end

[dummySort,indicesSorted]=sort([allPairs.diff]);
allPairs = allPairs(indicesSorted);
%%
% Initiate Union-Find data structure
global unionFind;
unionFind(sizeY,sizeX).parent.x = sizeX;
unionFind(sizeY,sizeX).parent.y = sizeY;
unionFind(sizeY,sizeX).size = 1;
unionFind(sizeY,sizeX).avgDx = dxs(sizeY,sizeX);
unionFind(sizeY,sizeX).avgDy = dys(sizeY,sizeX);
for y = 1 : sizeY 
    for x = 1 : sizeX
        unionFind(y,x).parent.x = x;
        unionFind(y,x).parent.y = y;
        unionFind(y,x).size = 1;
        unionFind(y,x).avgDx = dxs(y,x);        
        unionFind(y,x).avgDy = dys(y,x);        
    end
end
%%
% Iterate over pairs and merge pairs if they hold the predicate
for i = 1 : size(allPairs,2)
    curPair = allPairs(i);
    [region1] = findRegion(curPair.y1,curPair.x1);
    [region2] = findRegion(curPair.y2,curPair.x2);
    if (region1.parent.x == region2.parent.x && region1.parent.y == region2.parent.y)
        continue;
    end
    v1 = sqrt(region1.avgDx^2 + region1.avgDy^2);
    v2 = sqrt(region2.avgDx^2 + region2.avgDy^2);
    dist = sqrt((region1.avgDx - region2.avgDx)^2 + (region1.avgDy - region2.avgDy)^2)/max(v1,v2);
    if (mergePredicate(dist,region1,region2,params))
        %[unionFind]  = merge(region1,region2,unionFind);
        merge(region1,region2);
    end
end

% arrange output image
for y = 1 : sizeY 
    for x = 1 : sizeX
        [region]  = findRegion(y,x);
        outImgDx(y,x) = region.avgDx;
        outImgDy(y,x) = region.avgDy;
    end
end
% figure; imagesc(outImgDx); caxis([-10 10]); title('RMS: output image Dx');
% figure; imagesc(outImgDy); caxis([-10 10]); title('RMS: output image Dy');

end

%%
%%
%function [region, unionFind] = findRegion(unionFind, y,x)
function [region] = findRegion(y,x)
global unionFind;
%[parent, unionFind] = findParent(unionFind, y,x);
[parent] = findParent(y,x);
region = unionFind(parent.y,parent.x);
end

%function [parent, unionFind] = findParent(unionFind, y,x)
function [parent] = findParent(y,x)
global unionFind;
curRegion = unionFind(y,x);
if ((curRegion.parent.y == y) && (curRegion.parent.x == x))
    parent = curRegion.parent;    
else
    %[parent, unionFind] = findParent(unionFind,curRegion.parent.y,curRegion.parent.x);
    [parent] = findParent(curRegion.parent.y,curRegion.parent.x);
    unionFind(y,x).parent = parent;
end
end

function[merge] = mergePredicate(diff,region1,region2,params)
b1 = b(region1,params);
b2 = b(region2,params);
if (diff < (b1 + b2))
    merge =  1;
else 
    merge = 0;
end
end

%function[unionFind] = merge(region1,region2,unionFind)
function [] = merge(region1,region2)
global unionFind;
if (region1.size < region2.size)
    tmp = region1;
    region1 = region2;
    region2 = tmp;
end
% merge region2 to region1
mergeSize = region1.size + region2.size;
region1.avgDx = (region1.avgDx * region1.size + region2.avgDx * region2.size) / mergeSize;
region1.avgDy = (region1.avgDy * region1.size + region2.avgDy * region2.size) / mergeSize;
region1.size = mergeSize;
unionFind(region1.parent.y,region1.parent.x) = region1;% update region1 with the total region
unionFind(region2.parent.y,region2.parent.x).parent = region1.parent; % region2 parent points to region1
end

function [score] = b(region,params)
% global Q;
global contLP;
global logTable;

% score = sqrt(contLP + logTable(region.size + 1))/(2*Q*region.size);
score = (contLP * logTable(region.size + 1) * params.Q);
end