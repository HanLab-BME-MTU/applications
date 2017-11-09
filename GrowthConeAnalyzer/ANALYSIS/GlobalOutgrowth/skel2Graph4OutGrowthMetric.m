function [adjMat,edges,tipVertices,tipXYCoords,vertMat,edgePixels,error] = skel2Graph4OutGrowthMetric( skelIn )
%
% INPUT: some form of skeleton 
% OUTPUT: 
% adjMat: the matrix to be put into the graphshortestpath function 
% edges: nx2 array giving the vertices at the end of each possible edge
% 
% tipVertices: tip specific vertices for input into the graphshortest path function   
% vertMat: a label matrix where the vertices are labeled by ID 
% edgePixels: the actual pixel values between edges (ordered along the path so can make accurate dist measurements)
% error: scalar 1 or 0 : 1 if something was funny with the morphology
%% CHECK INPUT 
% 
if nargin< 1|| isempty(skelIn) || ~islogical(skelIn) || ndims(skelIn) ~=2
    error(':0 The first input must be a 2D binary matrix containing a skeleton: Try Again :0') 
end 
%%
% Clean up the skeleton 
% CC = bwconncomp(skelIn);
% csize = cellfun(@(c) numel(c), CC.PixelIdxList); 
% % get rid of anything that is not the largest CC 
% nsmall = sum(csize<max(csize)); 
% CC.NumObjects = CC.NumObjects - nsmall; 
% CC.PixelIdxList(csize<max(csize)) = []; 
% skelIn = labelmatrix(CC)>0; 

%endpoints = bwmorph(skelIn,'endpoints'); % output is mask of endpoints
%Note bwmoroph endpoints actual givens some incorrect points. 
% therefore directly find nearest neighbor of 1 here.
%note having a little of trouble with this nov 26th so returned to fine
%enpoints via bwmorph
retestSkel = 1; 
while retestSkel == 1 
nn = padarrayXT(double(skelIn~=0), [1 1]);
sumKernel = [1 1 1];
nn = conv2(sumKernel, sumKernel', nn, 'valid');
nn1 = (nn-1) .* (skelIn~=0);
edgeMat = nn1 == 2; % get pixels corresponding to edges 
vertMat = nn1 > 2| nn1==1; 
% maybe just get the endpoints from the edgeMat

nnEdge = padarrayXT(double(edgeMat~=0), [1 1]);
sumKernel = [1 1 1];
nnEdge = conv2(sumKernel, sumKernel', nnEdge, 'valid');
nn1Edge = (nnEdge-1) .* (edgeMat~=0);
endpointsEdge = find(nn1Edge==1); % will have overkill (ie calc each edge by def twice) 
% however could just find unique as a quick fix. 

nendpoints = length(endpointsEdge); 
 
% get labels 
%[edgeMat,nEdges] = bwlabeln(edgeMat); 
[vertMat, nVerts] = bwlabeln(vertMat); 

[ny,nx] = size(skelIn); 

for j = 1:nVerts
    
    %Get the index of the point(s)
    currInd = find(vertMat == j);
    
    %Convert this to matrix coord
    [currM,currN] = ind2sub(size(skelIn),currInd);
    
    %Average/store these coord
    vertices(j,:) = [mean(currM),mean(currN)];        


end

% get length of edges 

%edgePixels =  CC.PixelIdxList; 
% NOTE 2014_02_04 this is not completely correct need to convert to x,y coords and calc dist  
mkPlot2 =0;
if mkPlot2 == 1
    figure; 
    imshow(skelIn,[]) 
    hold on 
end 
% problem is that the data is not in order!! 
% need to check if see only one CC! else going to be double counting your
% paths. 
 CCTest = bwconncomp(edgeMat>0); %
% % get rid of singletons % bad idea 20140214
% singleOut = cellfun(@(x) length(x)==1,CCTest.PixelIdxList); 
% CCTest.PixelIdxList(singleOut) = []; 
% CCTest.NumObjects = CCTest.NumObjects- sum(singleOut); 

%Note not a good idea to remove all singletons as some are spanning the
%tree solution fixed getEndpoints to recognized singletons 


% if CCTest.NumObjects ==1 
%     % remove one of the endpoints 
%     % just take the first endpoint as this is corresponding to two sides. 
%     endpoints =  endpoints(1); 
%     nendpoints = 1; 
% end 
% get the endpoints of the edges 
EPsEdges  = cellfun(@(x) getEndpoints(x,[ny,nx],1),CCTest.PixelIdxList,'uniformoutput',0); % flag for singletons
% just take one (all need to make sure have in order) 
EPsEdges = cellfun(@(x) x(1,:),EPsEdges,'uniformoutput',0); 

nendpoints = numel(EPsEdges); 
CC.PixelIdxList = cell(1,nendpoints);
 for j = 1:nendpoints
%     
%     %Get the index of the point(s)
     %currInd = find(endpointMat == j);    
%     %Convert this to matrix coord (row is y)
% get coords 

% 
    pixIdxBack = nan(50,1); % overinitialize to make happy
%    
    transform = bwdistgeodesic(edgeMat,EPsEdges{j}(:,1),EPsEdges{j}(:,2)); % input in col(x),row(y) (xy coords)
%  
    iPix = 0; 
   while length(find(transform==iPix)) == 1
      pixIdxBack(iPix+1) = find(transform==iPix); % start at the endpoint
      iPix = iPix +1;  
   end 
%   
    pixIdxBack = pixIdxBack(~isnan(pixIdxBack)); 
 % make own pixelIdxList 
 CC.PixelIdxList{j} = pixIdxBack;
 
end 


edgePixels = CC.PixelIdxList; 

edgeLen = cellfun(@(x) calculateDistance(x,[ny,nx],mkPlot2) , CC.PixelIdxList)';


% CC.PixelIdxList(edgeLen==1) = []  ; 
% 
% CC.NumObjects = CC.NumObjects-sum(edgeLen==1); 

%edgeLen= cellfun(@(x) length(x),CC.PixelIdxList)'; % old incorrect way to



% sometimes the labels aren't corresponding between bwlabeln an the
% bwconncomp: so make own labelMat here to make sure edgeLen corresponds 
edgeMatLabel = zeros(size(edgeMat));
nEdges = numel(CC.PixelIdxList);
for iEdge = 1:numel(CC.PixelIdxList) 
edgeMatLabel(CC.PixelIdxList{iEdge}) = iEdge; 

end 

%badEdgeFlag = 0; 
for iEdge = 1:nEdges    
     %Find vertices which this edge connects
     % dilate the edge and find overlapping labels of vertices 
     x =unique(vertMat(imdilate(edgeMatLabel == iEdge,strel('disk',2)) & vertMat)); 
 if length(x) == 1 
      x = [x,x]; 
 end 

% if doesn't span two vertices likely too short-  just noise 
% if length(x) ~=2
%    % delete the edge 
%    badEdgeFlag = badEdgeFlag +1; 
%    edges(iEdge,:) = [NaN,NaN]; 
% else 

   
   

 edges(iEdge,:) = x; 
%end 

end % iEdge 
% if badEdgeFlag >0
%   edges(isnan(edges(:,1)),:) = [];
% 
% nEdges = nEdges - badEdgeFlag; 
% end 


% you need to test to make sure that the branch ID and the endpoint ID are
% not identical (might have small endpoints of 1 pixel near a branch) 

endpoints = nn1==1;

% find the vertices that are endpoints 
tipVertices =  vertMat.*endpoints; 
tipVertices = tipVertices(:); 
tipVertices = unique(tipVertices(tipVertices~=0));

% check to make sure tipVertex and branchvertex doesn't have the same label
% 
branchPoints = nn1 > 2; 
branchVertices = vertMat.*branchPoints; 
branchVertices = branchVertices(:); 
branchVertices = unique(branchVertices(branchVertices~=0)); 

allVert = [tipVertices;branchVertices]; 
 noRepeat = unique(allVert); 
 if length(noRepeat) == length(allVert)
     retestSkel = 0; 
 else 
 % fix the vertices so that remove the small endpoint and branch and
 % recalculate. 
 %toremove = setdiff(noRepeat,allVert); 
 % find the repeat 
 num = arrayfun(@(x) sum(allVert==noRepeat(x)),1:length(noRepeat)); 
 tipVertices = vertMat.*endpoints; 
 vertToRemove  = noRepeat(num>1);
 for iVert = 1:length(vertToRemove) 
 vertMat(vertMat == vertToRemove(iVert))=0;
  skelIn(tipVertices==vertToRemove(iVert)) = 0; 
 end 
 
% get rid of it from the skeleton 
 % recalc edges 
 clear edges edgeLen 
 end 
end % while 


% create the adjacency matrix (need to create bi-directional paths between nodes) 

adjMat =  sparse(vertcat(edges(:,1),edges(:,2)),vertcat(edges(:,2),edges(:,1)),repmat(edgeLen,[2 1]),nVerts,nVerts,2*nEdges); 







tipXYCoords = zeros(length(tipVertices),2); 
% check for end points and  branch points that have the same labels. 


% get the coordinates that correspond to these endpoints 
for i = 1:length(tipVertices) 
    idx= find(vertMat==tipVertices(i)); 
    
        
    [yTip,xTip] = ind2sub(size(skelIn),idx); 
    
    tipXYCoords(i,:) = [xTip, yTip]; 
end 

error = 0; 
% get the coordinates of these vertices 



% branchpoints = bwmorph(skelIn,'branchpoints'); % output is mask of branchpoints


% get x,y coords of end points
% [endpointMat, nendpoints] = bwlabeln(endpoints);  
% [branchpointMat,nbranchpoints] = bwlabeln(branchpoints); 
% 
% 
% 
% 
% branchPointMask = zeros(size(skelIn));
% % eventually incorporate into hunters skel2graph function for 3D 
% verticesEP = zeros(nendpoints,2); 
% verticesBP = zeros(nbranchpoints,2); 
% edgePathCoord = cell(nendpoints,1); 

% for j = 1:nendpoints
%     
%     %Get the index of the point(s)
%     currInd = find(endpointMat == j);    
%     %Convert this to matrix coord (row is y)
%     [currM,currN] = ind2sub(size(skelIn),currInd);
%     %Average/store these coord
%     verticesEP(j,:) = [floor(mean(currM)),floor(mean(currN))];   
% 
%     pixIdxBack = nan(50,1); % overinitialize to make happy
%    
%     transform = bwdistgeodesic(skelIn,verticesEP(j,2),verticesEP(j,1));
%  
%     iPix = 1; 
%   while length(find(transform==iPix)) == 1
%      pixIdxBack(iPix) = find(transform==iPix); % 
%      iPix = iPix +1;  
%   end 
%   
%    pixIdxBack = pixIdxBack(~isnan(pixIdxBack)); 
%     
%    if ~isempty(pixIdxBack)
%        possibleBranchPt{j} =  [find(transform == iPix) ; pixIdxBack(end)];
%    else
%        possibleBranchPt{j} = find(transform ==iPix);
%    end
%    branchPointMask(possibleBranchPt{j}) = 1;
%   % include last coord in possible branch points
%   
%   
%    
%     
%              
%     % get coords of the these pixels 
%     [yBack,xBack]= ind2sub(size(skelIn),pixIdxBack);
%     if ~isempty(yBack)
%     verticesBP(j,1) = yBack(end); 
%     verticesBP(j,2) = xBack(end);
%     edgePathCoord{j} = [yBack,xBack];
%     
%     else 
%         verticesBP(j,1) = verticesEP(j,1); 
%         verticesBP(j,2) = verticesEP(j,2); 
%         edgePathCoord{j} = verticesEP(j,:); 
%     end 
%    
%    
%     
%     
% end

% % test to see if the branchpoints are the same 
% for j = 1:nbranchpoints
%     
%     %Get the index of the point(s)
%     currInd = find(branchpointMat == j);
%      
%     %Convert this to matrix coord (row is y)
%     [currM,currN] = ind2sub(size(skelIn),currInd);
%     
%     %Average/store these coord
%     verticesBP2(j,:) = [floor(mean(currM)),floor(mean(currN))];      
% end


end

