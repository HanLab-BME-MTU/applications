function [neuriteLength,longPathLinInd,EPLongPath,error] =  GCAfindVeilStemLongestPath(veilStemMaskC,neuriteEntranceCLinIdx)
%% GCAfindVeilStemLongestPath:
% From a veilStemMask (a mask of the neurite/growth cone minus
% thin protrusions like filopodia)- this function skeletonizes the veil/stem mask
% and finds the longest path through the this skeleton starting from the
% point of the veil/stems entrance into the image.
% Here we are assuming that fluctuations in thin bundles of filamentous actin
% do not contribute to the overall neurite length.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%
% veilStemMaskC: (REQUIRED) RxC logical array
%        where R is the height (ny) and C is the width
%       (nx) of the input image (Output of Step III of GCA Segmentation)
%
% neuriteEntranceCLinIdx (REQUIRED) 1x1 double array
%        marking the linear index of the neurite entrance (Output of Step
%        III of GCA Segmentation)
%
% img : (OPTIONAL)
%
% OUTPUT:
% neuriteLength: (1x1) double
%   the neurite length measurement in um
%
% longPathLinInd: rx1 double
%   Where r is the number of pixels in the longest path
%   Output in linear index format
%
% EPLongPath: 1x2 double
%   r,c coords of endpoint of the longest path - used as input in 
%% PARAMS
 
%% CHECK INPUT

%% Initiate
error = false;
[ny,nx] = size(veilStemMaskC); 

%%   % I think that the 2D skeletonization method is more appropriate here.
%thinnedBodyAll = bwmorph(veilStemMaskC,'skel','inf');
%% this is internal for didatic purposes to show the skeletonization once

makeSmallMovie =0;
if makeSmallMovie == 1
    thinnedBodyAll = veilStemMaskC;
    same = 0 ;
    count = 1;
    roiYX = bwboundaries(veilStemMaskC);
    
    while same ==0
        thinnedBodyAll1 = bwmorph(thinnedBodyAll,'skel',count);
        imagesc(bwdist(~veilStemMaskC));
        imshow(thinnedBodyAll,[]);
        hold on
        cellfun(@(x) plot(x(:,2),x(:,1),'color','w'),roiYX);
        same = isequal(thinnedBodyAll1,thinnedBodyAll);
        thinnedBodyAll = thinnedBodyAll1;
        
        
        saveas(gcf,[num2str(count,'%03d') '.tif']);
        count = count+1;
        close gcf
    end
else
    thinnedBodyAll = bwmorph(veilStemMaskC,'skel','inf');
end

%% CREATE GRAPH FROM SKELETON
[adjMat,edges , tipVertices,tipXYCoords,vertMat,edgePixels]   = skel2Graph4OutGrowthMetric(thinnedBodyAll);

%% FIND THE SKELETON TIP THAT IS CLOSEST TO THE NEURITE COORD %%

[yEnter,xEnter] = ind2sub(size(thinnedBodyAll),neuriteEntranceCLinIdx);
%    end
nTips = length(tipXYCoords(:,1));

%get distances of each tip to the enter coord
distFromEnt = arrayfun(@(i) sqrt((xEnter-tipXYCoords(i,1))^2 + (yEnter-tipXYCoords(i,2))^2),1:nTips);

% find tip with teh min distance: this will be start point for the graph
% shortest path algorithm
enterVertex = tipVertices(distFromEnt==min(distFromEnt));

% the rest of the skeleton tips will be the other input into the graphshortestpath algorithm
tipVertices(distFromEnt==min(distFromEnt)) = []; % take out the enter index

if ~isempty(tipVertices) % % implement a small sanity check - if tipVerices found continue
    
    %% FIND THE LONGEST PATH OF THE SKELETON FROM THE NEURITE ENTRANCE TO THE END OF THE GROWTH CONE
    % This function will find the shortest path between the neurite enter vertex
    % and each of the tips of the neurite skeleton on the growth cone given the
    % distance edge weights
    % from there one can choose the longest of these paths as a good approximation of the neurite
    % length metric for a given frame
    for iTip = 1:length(tipVertices)
        [distPath(iTip), path{iTip} ] = graphshortestpath(adjMat,enterVertex,tipVertices(iTip),'Directed',false);
    end
    
    % Find max distance path
    idxMax =  find(distPath==max(distPath));
    idxMax  = idxMax(1); % just choose the first in case there is a tie.
    longPath = path{idxMax};
    %% MAKE FINAL LONG PATH MASKS AND RECALCULATE THE FINAL DISTANCE OF THE LONG PATH
    % NOTE: in the original implementation the length to the vertices is not included in the final
    %       distPath calc therefore it is slightly more correct to remake the path mask
    %       with the vertices and re-calculate the distance
    %       -20140425 should be fixed check this maria and give final ok before
    %       release
    
    edgesAll = [vertcat(edges(:,1), edges(:,2)), vertcat(edges(:,2),edges(:,1))];
    idxEdge = zeros(length(longPath)-1,1);
    % plot the pixels that correspond to the path lenghs
    for iPath = 2:length(longPath)
        startNode = longPath(iPath-1);
        endNode = longPath(iPath);
        %test(i-1,:) = [startNode,endNode];
        idxLog = arrayfun(@(x) isequal(edgesAll(x,:),[startNode,endNode]),1:length(edgesAll(:,1)));
        idxEdge(iPath-1) = find(idxLog==1);
    end
    edgePixels = [edgePixels';edgePixels'] ;
    
    
    
    % get edge Pixels
    pathMask = zeros(size(thinnedBodyAll));
    pathMask(vertcat(edgePixels{idxEdge})) = 1;
    
    
    idxVert = arrayfun(@(x) find(vertMat==x),longPath,'uniformoutput',0);
    % add the path vertices to the path mask
    
    pathMask(vertcat(idxVert{:})) = 1;
    
    thinMask = 1;
    while thinMask == 1
        pathMask   = bwmorph(pathMask,'thin','inf');
        pathMask = bwmorph(pathMask,'spur');
        pathMask(vertcat(idxVert{[1,end]})) = 1; % make sure to add back the vertices
        
        % keep testing the path mask for junctions
        nn = padarrayXT(double(pathMask~=0), [1 1]);
        sumKernel = [1 1 1];
        nn = conv2(sumKernel, sumKernel', nn, 'valid');
        nn1 = (nn-1) .* (pathMask~=0);
        
        test= nn1 > 2;
        % no more problem points you are done
        if sum(test(:)) ==0
            thinMask = 0;
        end
    end
    
    finalPathPix= find(pathMask==1);
    EPLongPath = getEndpoints(finalPathPix,[ny,nx],1);
    [nEndpoints,~]= size(EPLongPath);
    dist = arrayfun(@(i) sqrt((xEnter-EPLongPath(i,1))^2 + (yEnter-EPLongPath(i,2))^2),1:nEndpoints);
    
    
    EPLongPath = EPLongPath(dist==max(dist),:); % take the end farthest away from the entrance
    
    % Reorder the pixels
    longPathLinInd = nan(length(finalPathPix),1); % overinitialize to make happy
    transform = bwdistgeodesic(pathMask,EPLongPath(:,1),EPLongPath(:,2)); % input in col(x),row(y) (xy coords)
    iPix = 0;
    while length(find(transform==iPix)) == 1
        longPathLinInd(iPix+1) = find(transform==iPix); % start at the endpoint
        iPix = iPix +1;
    end
    longPathLinInd = longPathLinInd(~isnan(longPathLinInd));
 
    
    % Recalculate the distance
    [ neuriteLength]  = calculateDistance(longPathLinInd,[ny,nx],0);
    
    
    
    
    
else
    display(['Error: No Skeleton Vertices Found']);
    error = true ;
    
end % ~isempty(tipVertices




end

