function [verticesEP,verticesBP, edgePathCoord,branchPointMask] = skel2graph2D( skelIn )
%
% 2D version of the hunters function skel2graph3D need to eventually
% incorporate!!!!!
% OUTPUT:
% verticesEP: m x 2 matrix (m being the number of enpoints) containing the
% coordinates cooresponding to the endpoints of the skeleton
%
% verticesBP: mx2 matrix of "Branch" points of skel:  for now just the idx of the pixel before
% dist transform is non-unique > 1 value (maybe change in future)
% Note to self: look up officially def of branch ...
%
% edgePathCood: m cell array containing the yx coordinates of all paths
% between two vertices connected by an "edge"
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
nn = padarrayXT(double(skelIn~=0), [1 1]);
sumKernel = [1 1 1];
nn = conv2(sumKernel, sumKernel', nn, 'valid');
nn1 = (nn-1) .* (skelIn~=0);
edgeMat = nn1 == 2; % get pixels corresponding to edges
vertMat = nn1 > 2| nn1==1;

% get labels
[edgeMat,nEdges] = bwlabeln(edgeMat);
[vertMat, nVerts] = bwlabeln(vertMat);

for j = 1:nVerts
    
    %Get the index of the point(s)
    currInd = find(vertMat == j);
    
    %Convert this to matrix coord
    [currM,currN] = ind2sub(size(skelIn),currInd);
    
    %Average/store these coord
    vertices(j,:) = [mean(currM),mean(currN)];
    
    
end


% for iEdge = 1:nEdges
%      %Find vertices which this edge connects
%      % dilate the edge and find overlapping labels of vertices
%  edges{iEdge}=  unique(vertMat(imdilate(edgeMat == iEdge,strel('disk',2)) & vertMat));
% end

% if mkPlot == 1
%      spy(skelIn,'r')
%      hold on
%      % plot the vertice with it's label
%      for iEdge = 1:nEdges
%      [y,x] = ind2sub(size(skelIn),find(vertMat==iEdge));
%      scatter(y(1),x(1), 10,'y')
%
%
%      end
endpoints = nn1==1;

branchpoints = bwmorph(skelIn,'branchpoints'); % output is mask of branchpoints


% get x,y coords of end points
[endpointMat, nendpoints] = bwlabeln(endpoints);
[branchpointMat,nbranchpoints] = bwlabeln(branchpoints);




branchPointMask = zeros(size(skelIn));
% eventually incorporate into hunters skel2graph function for 3D
verticesEP = zeros(nendpoints,2);
verticesBP = zeros(nbranchpoints,2);
edgePathCoord = cell(nendpoints,1);

for j = 1:nendpoints
    
    %Get the index of the point(s)
    currInd = find(endpointMat == j);
    %Convert this to matrix coord (row is y)
    [currM,currN] = ind2sub(size(skelIn),currInd);
    %Average/store these coord
    verticesEP(j,:) = [floor(mean(currM)),floor(mean(currN))];
    
    pixIdxBack = nan(50,1); % overinitialize to make happy
    
    transform = bwdistgeodesic(skelIn,verticesEP(j,2),verticesEP(j,1));
    % test
    %     figure;
    %     imagesc(transform)
    %     hold on
    %     scatter(verticesEP(j,2),verticesEP(j,1),10,'r','filled')
    
    iPix = 1;
    while length(find(transform==iPix)) == 1
        pixIdxBack(iPix) = find(transform==iPix); %
        iPix = iPix +1;
    end
    
    pixIdxBack = pixIdxBack(~isnan(pixIdxBack));
    
    if (~isempty(pixIdxBack) && length(pixIdxBack)>1)
        possibleBranchPt{j} = pixIdxBack(end); %branchPointMask(possibleBranchPt{j}) = 1
        %possibleBranchPt{j} =  [find(transform == iPix) ; pixIdxBack(end)];
        branchPointMask(vertcat(possibleBranchPt{j})) = 1;
    else
        possibleBranchPt{j} = [];
        % possibleBranchPt{j} = find(transform ==iPix);
    end
    
    
    %branchPointMask(vertcat(possibleBranchPt{:})) = 1;
    %branchPointMask(possibleBranchPt{j}) = 1;
    % include last coord in possible branch points
    % spy(branchPointMask,'w');
    
    
    
    
    % get coords of the these pixels
    [yBack,xBack]= ind2sub(size(skelIn),pixIdxBack);
    if ~isempty(yBack)
        verticesBP(j,1) = yBack(end);
        verticesBP(j,2) = xBack(end);
        edgePathCoord{j} = [yBack,xBack];
        
    else
        verticesBP(j,1) = verticesEP(j,1);
        verticesBP(j,2) = verticesEP(j,2);
        edgePathCoord{j} = verticesEP(j,:);
    end
    
    
    
    
end

% tentative : prune the junction mask
CCBs = bwconncomp(branchPointMask);
l = cellfun(@(x) length(x)>1,CCBs.PixelIdxList);

if sum(l)>0
    % get the first first pixel of each CC with more than one pixel
    toFix = find(l);
    for i = 1:length(toFix)
        pixDelete = CCBs.PixelIdxList{toFix(i)}(2:end);
        branchPointMask(pixDelete) = 0;
    end
    
end
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

