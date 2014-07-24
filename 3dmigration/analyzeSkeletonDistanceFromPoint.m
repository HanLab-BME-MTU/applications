function [edgeDistances,closestPt,dToClosest] = analyzeSkeletonDistanceFromPoint(vertices,edges,edgePaths,queryPt)
%ANALYZESKELETONDISTANCEFROMPOINT finds distance ALONG THE SKELETON from the input point at each point on the skeleton
%
% Query pt should be in image coordinates (y,x,z) ! 
%
% TEMP - ADD MORE DOCUMENTATION!!!
% TEMP - this function doesn't handle loops well!! Either fix it, or open
% the loops before it's called! As long as the query point is close to the
% center of the cell, the edge that is excluded is usually the
% false-positive anyways, so for now the error introduced is minor...
% 
%   For per-vertex distance from tip to skeleton body, see
%   analyzeSkeletonTipPaths.m
%
% Hunter Elliott
% 2/2012
%

showPlots = false;%For testing/debugging

nEdges = numel(edgePaths);
nVert = size(vertices,1);
nPtsPerEdge = cellfun(@(x)(size(x,1)),edgePaths);
edgeDirs = getEdgePathDir(vertices,edges,edgePaths);


% ---- Find closest point on the skeleton to query pt ---- %


%Combine all the skeleton coordinates so we can create the KDTree only once
allEdgePts = vertcat(edgePaths{:});
indArray = arrayfun(@(x)(repmat(x,nPtsPerEdge(x),1)),1:nEdges,'Unif',false);%Allows us to determine which edge the closest point came from
indArray = vertcat(indArray{:});

%Now find the closest point and its distance
[iQueryEdgePt,dToClosest] = KDTreeClosestPoint(allEdgePts,queryPt);

%Get the edge to which it belongs, which we call the query edge
closestPt(1) = indArray(iQueryEdgePt);

%..and the index of the point on the query edge
tmp = false(size(indArray));
tmp(iQueryEdgePt) = true;
closestPt(2) = find(tmp(indArray == closestPt(1)));

%..and the vertices corresponding to this edge
queryEdgeVert = edges(closestPt(1),:);


% ---- Calc distance along skeleton to this point from every point ---%

edgeDistances = cell(nEdges,1);
edgePtSpacing = cellfun(@(x)(sqrt(sum(diff(x,1,1) .^2,2))),edgePaths,'Unif',0);%Spacing between edge points
edgeLen = cellfun(@sum,edgePtSpacing);%And the length of each edge
%Due to the image discretization and skeleton-graphization method, some
%edges are a single point. Set the length to 1 so that these will still
%provide connectivity in the adjacency matrix.
edgeLen(edgeLen==0) = 1;



%Create adjacency matrix with edge lengths as weights
adjMat = sparse(vertcat(edges(:,1),edges(:,2)),vertcat(edges(:,2),edges(:,1)),repmat(edgeLen,[2 1]),nVert,nVert,2*nEdges);

%The edge containing the query point is a special case - do it first;
tmp1 = cumsum(edgePtSpacing{closestPt(1)}(closestPt(2)-1:-1:1));
tmp2 = cumsum(edgePtSpacing{closestPt(1)}(closestPt(2):end));
edgeDistances{closestPt(1)} = vertcat(tmp1(end:-1:1),0,tmp2);

%TEMP - This is not the fastest way to do this (some redundancy in edge
%checking)... BUt plenty fast enough
for j = 1:nVert
    
    %Get the edges for this vertex
    iCurrEdges = any(bsxfun(@eq,edges,j),2);    

    if any(cellfun(@isempty,edgeDistances(iCurrEdges)))
       
        %Find the minimum distance to each of the vertices on the query
        %edge
        currPathVert = cell(2,1);
        currPL = nan(2,1);
        for k = 1:2
            [currPL(k) currPathVert{k}] = graphshortestpath(adjMat,j,queryEdgeVert(k),'Directed',false);
        end
        %First make sure ther is a valid path - image edge effects aren't
        %handled well in the current skeleton-graphization so there are
        %occasionally isolated fragments.... TEMP!
        if any(isfinite(currPL))
            [~,iCurrCloseEdgeVert] = min(currPL);
            currPathVert = currPathVert{iCurrCloseEdgeVert};
            nPtsCurr = numel(currPathVert);
            currPathEdges = nan(nPtsCurr-1,1);
            %Now find the edges which compose this path... Convert this to a
            %bsxfun call?
            for k = 1:nPtsCurr-1
                %TEMP - If this edge has multiplicity > 1 take the shortest edge,
                %since we are finding shortest-path anyways. Eventually, we
                %should open loops before even calling this function.
                tmpEdges = find(any(edges == currPathVert(k),2) & any(edges == currPathVert(k+1),2));
                if numel(tmpEdges) > 1
                    [~,iShortest] = min(edgeLen(tmpEdges));
                    tmpEdges = tmpEdges(iShortest);       
                end
                currPathEdges(k) = tmpEdges;

            end 

            %First get the distance from the end of this path to the query
            %point - we will need ot add this to each edge path length
            if currPathVert(end) == edges(closestPt(1),1)
                hitsStart = 1;
            else
                hitsStart = -1;
            end
            %Depending on the direction of the edge path and the vertex we hit,
            %we sum either the first or last section of the query edge as our
            %offset
            if hitsStart*edgeDirs(closestPt(1)) > 0 
                geOffset = sum(edgePtSpacing{closestPt(1)}(1:min(closestPt(2),end)));
            else
                geOffset = sum(edgePtSpacing{closestPt(1)}(min(closestPt(2),end):end));
            end


            %Now go through each edge on this path and calculate it's per-point
            %distance
            for k = 1:numel(currPathVert)-1

                %Get the distance from the end of this edge to the query pt            
                dOffset = sum(edgeLen(currPathEdges(k+1:end)));%Add the distance from the next vertex to the first vertex on the query edge...            
                dOffset = dOffset + geOffset;%And the distance from that query vertex to the query point
                %Depending on the direction of the path relative to the edge
                %path and edge vertices, we need to reverse the cumulative
                %sum...
                edgeOrdEqPathOrd = all(currPathVert(k:k+1) == edges(currPathEdges(k),:));            
                if (edgeOrdEqPathOrd && edgeDirs(currPathEdges(k))==1) || (~edgeOrdEqPathOrd && edgeDirs(currPathEdges(k))==-1)
                    edgeDistances{currPathEdges(k)} = vertcat(0,cumsum(edgePtSpacing{currPathEdges(k)}(end:-1:1)));
                    edgeDistances{currPathEdges(k)} = edgeDistances{currPathEdges(k)}(end:-1:1);
                else
                    edgeDistances{currPathEdges(k)} = vertcat(0,cumsum(edgePtSpacing{currPathEdges(k)}));
                end
                %Add the offset
                edgeDistances{currPathEdges(k)} = edgeDistances{currPathEdges(k)} + dOffset;            
            end
        end
    end        
    
end
    
if showPlots
    figure    
    hold on
    plot3(queryPt(2),queryPt(1),queryPt(3),'rx','MarkerSize',15);    
    plot3(allEdgePts(iQueryEdgePt,2),allEdgePts(iQueryEdgePt,1),allEdgePts(iQueryEdgePt,3),'gx','MarkerSize',15)
    plot3(edgePaths{closestPt(1)}(closestPt(2),2),edgePaths{closestPt(1)}(closestPt(2),1),edgePaths{closestPt(1)}(closestPt(2),3),'go','MarkerSize',15);
    legend('Query Pt','Closest','Closest By Skel')
    
    plotSkel(vertices,edgePaths);
    %plotDirection3(edgePaths);
    iHasDist = find(~cellfun(@isempty,edgeDistances));
    arrayfun(@(x)(scatter3(edgePaths{x}(:,2),edgePaths{x}(:,1),edgePaths{x}(:,3),edgeDistances{x}*5+1e-3)),iHasDist)
end
    
    
