function [normalsCRotated,smoothedEdgeC,normalsC] = gcaReorientVeilStemNormalsTowardsOutgrowth(leadProtrusionPt,LPIndices,normalsC,smoothedEdgeC,dims,idxEnter )
% gcaReorientVeilStemNormalsTowardOutgrowth
%
% leadProtrusionPt:
%
% normalsC:
%
% smoothedEdgeC:
%
% dims :
%
% OUTPUT:
% normalsCRotated
%% Initiate
nSamples = size(normalsC,1);
% initiate rotated normals mat
normalsCRotated = zeros(nSamples,2);

% find the closest point to the smoothed edge values as these will be our
% normals
%get distances of each tip to the enter coord
distToSmoothEdge = arrayfun(@(i) sqrt((leadProtrusionPt(1,1)-smoothedEdgeC(i,1))^2 + (leadProtrusionPt(1,2)-smoothedEdgeC(i,2))^2),1:nSamples);

% get the distances to the smoothed edge neurite entrance
[yEnter,xEnter] = ind2sub(dims,idxEnter);
distToSmoothEdgeNeuriteEnter = arrayfun(@(i) sqrt((xEnter-smoothedEdgeC(i,1))^2 + (yEnter-smoothedEdgeC(i,2))^2),1:nSamples);
% hold on
% scatter(smoothedEdgeC(:,1),smoothedEdgeC(:,2),10,'g','filled');
% hold on
% show the plot
% quiver(smoothedEdgeC(:,1),smoothedEdgeC(:,2),normalsC(:,1),normalsC(:,2),'g');
% scatter(leadProtrusionPt(:,1),leadProtrusionPt(:,2),'b');

% find the point on the smoothed edge that is closest to the lead
% protrusion
leadProtPtSmoothEdge= smoothedEdgeC(distToSmoothEdge == min(distToSmoothEdge),:);
neuriteEnterSmoothedEdge = smoothedEdgeC(distToSmoothEdgeNeuriteEnter ==min(distToSmoothEdgeNeuriteEnter),:);
%zeros(idxLead,)

% scatter(leadProtPtSmoothEdge(:,1),leadProtPtSmoothEdge(:,2),'y');

% create a boder mask to make sure that the length of the border = the the
% length of the smoothEdge coords which we will use for hte normal calcs
smoothedEdgeCIdxPix= sub2ind(dims,round(smoothedEdgeC(:,2)),round(smoothedEdgeC(:,1)));

borderMask = false(dims);
borderMask(smoothedEdgeCIdxPix) = true;
borderMaskAll = borderMask;

borderMask(round(leadProtPtSmoothEdge(1,2)),round(leadProtPtSmoothEdge(1,1))) = false;
borderMask(neuriteEnterSmoothedEdge(1,2),neuriteEnterSmoothedEdge(1,1)) = false;

%figure;

% the distMat doesn't always work well if the two sides not assymetric
% get the distance transform from the skeleto
%distMat  = bwdistgeodesic(borderMaskAll,round(leadProtPtSmoothEdge(1,1)),round(leadProtPtSmoothEdge(1,2)));
% distMat = distMat(1:max(distMat(:))-1);
%distMat(distMat==0 | distMat==max(distMat(:)))= NaN;



% maskD = false(dims);
% maskD(LPIndices) = 1;
% maskD = imdilate(maskD,strel('disk',1));
%
%
% distMat(maskD) = NaN;
% imshow(distMat>0,[]);
% hold on
% spy(maskD,'r',10);
% imagesc(distMat) ;
% get the two CCs created. (will be in wrong order)
CCs = bwconncomp(borderMask>0);


if CCs.NumObjects == 1
    % try to take the backbone mask and find the problem point
    backboneMask = false(dims);
    backboneMask(LPIndices) = true;
    backboneMask = bwmorph(backboneMask,'diag');
    toRemove = intersect(find(borderMaskAll==1),find(backboneMask==1));
    idxRemove = arrayfun(@(x)  find(smoothedEdgeCIdxPix == toRemove(x)),1:length(toRemove) );
    
    quiver(smoothedEdgeC(:,1),smoothedEdgeC(:,2),normalsC(:,1),normalsC(:,2),'g')
    normalsC(idxRemove,:) = [];
    smoothedEdgeC(idxRemove,:) = [];
    
    quiver(smoothedEdgeC(:,1),smoothedEdgeC(:,2),normalsC(:,1),normalsC(:,2),'r')
    smoothedEdgeCIdxPix(idxRemove) = [];
    normalsCRotated(idxRemove,:) = [];
    
    borderMask(backboneMask==1) = 0;
    borderMaskAll(backboneMask==1) = 0; 
    
    
    
    %distMat  = bwdistgeodesic(borderMask,round(leadProtPtSmoothEdge(1,1)),round(leadProtPtSmoothEdge(1,2)));
    %imagesc(distMat);
    %distMatBorder =  bwdist(borderMask);
    %imagesc(distMatBorder.*backboneMask);
    %LPIndices
    CCs = bwconncomp(borderMask>0);
    if CCs.NumObjects > 2 
       sizeSeg = cellfun(@(x) length(x),CCs.PixelIdxList);
       CCs.PixelIdxList(sizeSeg == min(sizeSeg)) = []; 
       CCs.NumObjects = CCs.NumObjects - 1;  
    end 
end

% make the distMat anyway just for testing the orientation relative to the
% protrusion
distMat  = bwdistgeodesic(borderMaskAll,round(leadProtPtSmoothEdge(1,1)),round(leadProtPtSmoothEdge(1,2)));

% match the indices of the smoothedEdge to the indices of the CCs

%smoothedEdgeCIdxPix = sub2ind(size(distMat),smoothedEdgeC(:,1),smoothedEdgeC(:,2));

% mark the pixIdxs to be rotated clockwise

[~,~,ib] = intersect(CCs.PixelIdxList{1},smoothedEdgeCIdxPix,'stable');

% mark the pixIdxs to be rotated counter

[~,~,ib2]= intersect(CCs.PixelIdxList{2},smoothedEdgeCIdxPix,'stable');
% figure
% make first mask
% mask1 = false(size(distMat));
% mask1(CCs.PixelIdxList{1}) = true;
% spy(mask1,'r');
% mask2 = false(size(distMat));
% hold on
% mask2(CCs.PixelIdxList{2}) = true;
% spy(mask2,'b')

%figure;
% figure out the directionality % if the normal should be rotatated + or  -
% clockwise or counterclockwise.
% if
% plot before
% quiver(smoothedEdgeC(:,1),smoothedEdgeC(:,2),normalsC(:,1),normalsC(:,2),'g');
%hold on

% test CCs dir
mask1 = false(dims);
mask1(CCs.PixelIdxList{1})= true;
CC1Dist= distMat.*mask1;
[y1,x1] = find(CC1Dist == 1);
% just take the first sometimes have more than one value 
x1 = x1(1); 
y1= y1(1); 

%find(smoothedEdgeCIdxPix,1);
idx1 = find(CC1Dist==1);
idx1 = idx1(1); 
toRot = [normalsC(smoothedEdgeCIdxPix==idx1,1),normalsC(smoothedEdgeCIdxPix==idx1,2)];
%make sure toRot is 
toRot = toRot(1,:); 

idx2 = find(CC1Dist==2);
[y2,x2] = find(CC1Dist == 2);
% just take the first (sometimes have more than one value) 
y2 = y2(1); 
x2 = x2(1); 

% scatter(x1,y1,'filled','r');
% scatter(x2,y2,'filled','g');
signs =[ -1,1,1,-1];


% direction toward leading prot
delty = (y1-y2) ;
deltx = (x1-x2);
vectTest1 = [deltx,delty];
% get the idx of

%quiver(smoothedEdgeC(ib(1),1),smoothedEdgeC(ib(1),2),vectTest1(1),vectTest1(2),'k')
rotMatrix = [0 -1 ; 1 0]; % 90 degree
% get the normals of the entrance to the
vectTest2 = toRot*rotMatrix;
%vectTest2 = idx2(1);
%quiver(smoothedEdgeC(idx,1),smoothedEdgeC(ib(1),2),vectTest2(1),vectTest2(2),'y')
test = dot(vectTest1,vectTest2);
% if the rotation if rotating away
if test<0
    signs = -(signs);
end

% plot after
rotMatrix = [0 signs(1) ; signs(2) 0]; % 90 degree
dirEdge = normalsC(ib,:)*rotMatrix;
normalsCRotated(ib,1:2) = normalsC(ib,:)*rotMatrix;
normalsCRotated(ib,3) = 1;

rotMatrix2 = [0 signs(3) ; signs(4) 0]; %
dirEdge2 = normalsC(ib2,:)*rotMatrix2;
normalsCRotated(ib2,1:2) = normalsC(ib2,:)*rotMatrix2;
normalsCRotated(ib2,3) = 2;

% sanity check
% side1 = find(normalsCRotated(:,3) == 1);
% side2 = find(normalsCRotated(:,3) == 2) ;
% quiver(smoothedEdgeC(side1,1),smoothedEdgeC(side1,2),normalsCRotated(side1,1),...
%     normalsCRotated(side1,2),'b');
% quiver(smoothedEdgeC(side2,1),smoothedEdgeC(side2,2),normalsCRotated(side2,1),...
%     normalsCRotated(side2,2),'r');

%quiver(smoothedEdgeC(ib(1),1),smoothedEdgeC(ib(1),2),vectTest1(1),vectTest1(2),'k');
%quiver(smoothedEdgeC(ib,1),smoothedEdgeC(ib,2),dirEdge(:,1),dirEdge(:,2),'r');


% hold on
% quiver(smoothedEdgeC(ib2,1),smoothedEdgeC(ib2,2),dirEdge2(:,1),dirEdge2(:,2),'b');



%idxStart = find(distMatInvert == max(distMatInvert(:))-1);
%[y,x] = ind2sub(size(distMatInvert),idxStart);

%[pts, values] = gradientDescent(distMatInvert,x(1),y(1));

%pts = gradientDescent(distMatInvert,skelPoint(1,1),skelPoint(1,2));
end
%% % reorder the CCS 

         
