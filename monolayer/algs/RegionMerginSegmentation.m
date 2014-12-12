function [outImg, unionFind] = RegionMerginSegmentation(inImg)
global P;
global Q;
global contLP;
global logTable;

inImg = double(inImg);
[sizeY,sizeX] = size(inImg);
outImg = zeros(sizeY,sizeX);

P = 0.00001;
Q = 0.001;
contLP = log(2.0/P);
logTable = log((1:(sizeY*sizeX) )+ 1);

% Initiate all pairs
numGroups = 2*(sizeY-2) * (sizeX-2);
allPairs(numGroups).x1 = 0;
allPairs(numGroups).y1 = 0;
allPairs(numGroups).x2 = 0;
allPairs(numGroups).y2 = 0;
allPairs(numGroups).diff = 0;
i = 1;
for y = 2 : sizeY - 1
    for x = 2 : sizeX - 1  
        % right
        allPairs(i).x1 = x;
        allPairs(i).y1 = y;
        allPairs(i).x2 = x+1;
        allPairs(i).y2 = y;        
        allPairs(i).diff = abs(inImg(y,x) - inImg(y,x+1));
        i = i + 1; % next group
        
        % down
        allPairs(i).x1 = x;
        allPairs(i).y1 = y;
        allPairs(i).x2 = x;
        allPairs(i).y2 = y+1;
        allPairs(i).diff = abs(inImg(y,x) - inImg(y+1,x));
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
unionFind(sizeY,sizeX).avgGL = inImg(sizeY,sizeX);
for y = 1 : sizeY 
    for x = 1 : sizeX
        unionFind(y,x).parent.x = x;
        unionFind(y,x).parent.y = y;
        unionFind(y,x).size = 1;
        unionFind(y,x).avgGL = inImg(y,x);        
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
    if (mergePredicate(abs(region1.avgGL - region2.avgGL),region1,region2))
        %[unionFind]  = merge(region1,region2,unionFind);
        merge(region1,region2);
    end
end

% arrange output image
for y = 1 : sizeY 
    for x = 1 : sizeX
        [region]  = findRegion(y,x);
        outImg(y,x) = round(region.avgGL);
    end
end
figure; colormap(gray); imagesc(outImg); title('RMS: output image');

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

function[merge] = mergePredicate(diff,region1,region2)
b1 = b(region1);
b2 = b(region2);
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
region1.avgGL = (region1.avgGL * region1.size + region2.avgGL * region2.size) / mergeSize;
region1.size = mergeSize;
unionFind(region1.parent.y,region1.parent.x) = region1;% update region1 with the total region
unionFind(region2.parent.y,region2.parent.x).parent = region1.parent; % region2 parent points to region1
end

function [score] = b(region)
global Q;
global contLP;
global logTable;

score = sqrt(contLP + logTable(region.size + 1))/(2*Q*region.size);
end