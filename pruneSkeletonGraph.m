function [vertices,edges,edgePaths,edgeLabels] = pruneSkeletonGraph(vertices,edges,edgePaths,mask,maskProp,varargin)
%PRUNESKELETONGRAPH prunes the input skeleton graph by geometric criteria
%
% [vertices,edges,edgePaths,labels] = pruneSkeletonGraph(vertices,edges,edgePaths,skelProps);
% [vertices,edges,edgePaths,labels] = pruneSkeletonGraph(...,'OptionName1',optionValue1,'OptionName2',optionValue2,...);
%
%   This function processes the input skeleton to remove, merge or split
%   the edges and vertices in the input skeleton based on several geometric
%   criteria, relying on various properties of the skeleton which must be
%   calculated using skelGraphProps.m Using these same geometric criteria,
%   the components of the skeleton are labelled as belonging to the "body"
%   or "branches" of the skeleton.
%
%   ***TEMP - THIS IS OUT OF DATE, AND NEEDS TO BE MORE DETAILED!!**** -
%   HLE
%   
%
% Input:
%
%   vertices,edges,edgePaths - The graph form of the skeleton, as returned
%   by skel2graph.m
%
%   skelProps - The skeleton properties, as returned by skelGraphProps.m
%
%   'OptionName',optionValue - A string with an option name followed by the
%   value for that option.
% 
%   Possible Option Names:
%
%       ('OptionName' -> possible values)
%
%       ('MaxRadius'->Positive scalar) The maximum radius (in voxels) of an
%       object to consider a branch. Skeleton branches resulting from areas
%       of the mask with radius larger than this will be pruned. Optional.
%       Default is 5.
%
%       ('MinAspectRatio'->positive scalar) Any portions of the skeleton
%       which correspond to areas of the mask having an aspect ratio
%       (length:radius) smaller than this value will be labelled as part of
%       the body of the skeleton, while the rest will be labelled branches
%       if they pass the other pruning criteria.
%       Optional. Default is 1.5
%
%       ('ShowPlots'->true/false) If true, figures will be displayed
%       showing the branches before and after pruning/labelling.
%       Optional. Default is false.
%   
% Output:
%
%    vertices,edges,edgePaths - The pruned skeleton graph. Same form as
%    input skeleton graph.
%
% Hunter Elliott
% 3/2011
%

%% ------------------- Parameters ---------------------- %%

nMultMax = 4;%Maximum length of edges with multiplicity > 1 to remove. See below for details.


%% ------------------------- Input -------------------------------%%

%Parse all the inputs
ip = inputParser;
ip.FunctionName = mfilename;

ip.addRequired('vertices',@(x)(ndims(x) == 2 && size(x,2)==3 && size(x,1)>1))
ip.addRequired('edges',@(x)(ndims(x) == 2 && size(x,2)==2 && size(x,1)>1))
ip.addRequired('edgePaths',@(x)(iscell(x) && numel(x) > 1));
ip.addRequired('mask',@(x)(islogical(x) && ndims(x) == 3));
ip.addRequired('maskProp',@(x)(isstruct(x) && numel(x) == 1));
ip.addParamValue('ShowPlots',false,(@(x)(numel(x)==1)));
ip.addParamValue('MaxRadius',5,(@(x)(numel(x) == 1 && x >= 1))); 
ip.addParamValue('MinAspectRatio',2,(@(x)(numel(x) == 1 && x >= 1))); 
ip.addParamValue('MinLength',5,(@(x)(numel(x) == 1 && x >= 1))); 
ip.addParamValue('ImageEdgeDist',3,(@(x)(numel(x) == 1 && x >= 1)));

ip.parse(vertices,edges,edgePaths,mask,maskProp,varargin{:});
p = ip.Results;


%TEMP - MORE CHECKING OF REQUIRED INPUTS???!!-HLE


%% -------------------------- Init ------------------------------%%

%Get radius-of-curvature thresholds from the maximum radius parameter

% gcThresh = 1/p.MaxRadius^2;%Maximum gaissian curvature value. This corresponds to
%                       %the gaussian curvature of a sphere of radius maxRad
%                          
% mcThresh = -1/p.MaxRadius; %Maximum mean curvature value. This corresponds to
%                       %the mean curvature of a sphere of radius maxRad

nVert = size(vertices,1);
nEdges = size(edges,1);
[M,N,P] = size(mask);

%% ----------------- Multiple-Edge Pruning ------------------- %%
%Removes short edges with multiplicity > 1 which are introduced by poor
%skeletonization/graphization

%First we deal with multiple-edges (two vertices connected by more than one
%edge) as these can be problematic later. 

%Calculate the multiplicity of each edge
edgeMultiplicity = nan(nEdges,1);
for j = 1:nEdges            
    edgeMultiplicity(j) = nnz(all(bsxfun(@eq,edges,edges(j,:)),2));    
end


iMult = find(edgeMultiplicity>1);
nMult = numel(iMult);
beenThere = false(nEdges,1);%Vector for keeping track of which edges we've dealt with
iMultPrune = nan(nMult,1);
for j = 1:nMult
    
    %Find which edge(s) share these vertices
    iSharing = iMult(any(bsxfun(@eq,edges(iMult,:),edges(iMult(j),:)),2));
    
    %Make sure we haven't dealt with this multiple-edge before
    if ~any(beenThere(iSharing))
        
        beenThere(iSharing) = true;
    
        %These are usually a result of the crappy skeletonization algortithm I'm
        %using, and are therefore usually very short edges with multiplicity =
        %2. These are the ones we want to remove. Longer edges or those with
        %higher multiplicity are a bad sign, so warn the user and leave them
        %alone.
        nPtsEach = cellfun(@(x)(size(x,1)),edgePaths(iSharing));
        if numel(iSharing) == 2 && all( nPtsEach <= nMultMax);

            %Merge the edges by replacing the first of the two edges with a
            %straight line
            startVert = vertices(edges(iSharing(1),1),:);
            endVert = vertices(edges(iSharing(1),2),:);
            edgePaths{iSharing(1)} = vertcat(startVert,endVert);        

            if p.ShowPlots   

                if ~exist('meFig','var')
                    meFig = plotSkel(vertices,edges,edgePaths); %#ok<NASGU>
                    title('Mutliple-Edge Pruning')
                end

                plot3(edgePaths{iSharing(2)}(:,2),edgePaths{iSharing(2)}(:,1),edgePaths{iSharing(1)}(:,3),'rx','MarkerSize',20)
                plot3(edgePaths{iSharing(1)}(:,2),edgePaths{iSharing(1)}(:,1),edgePaths{iSharing(1)}(:,3),'bo','MarkerSize',20)
            end

            %Mark the other edge for removal
            iMultPrune(j) = iSharing(2);
            


            %STOPPED HERE - HLE

        else
            warning('pruneSkel:edgeMultiplicity',...
                'Edge with multiplicity %d and max length %d was unable to be pruned!',numel(iSharing),max(nPtsEach));
        end           
    end
end

%Remove the multiple-edges which were marked for pruing
[vertices,edges,edgePaths,~,~] = deleteEdges(vertices,edges,edgePaths,iMultPrune(~isnan(iMultPrune)));

%Resolve any degree = 2 vertices. Some may have been created by the
%multi-edge pruning, others may have already existed in the skeleton. In
%either case, these are meaningless and problematic.
[vertices,edges,edgePaths,nEdges,nVerts] = removeDegree2Vert(vertices,edges,edgePaths,p);


%% -------- Minimum - Length Pruning ------- %%
%Removes terminal edges (tips) which are shorter than the minimum length

%TEMP - Should we loop through this to check length of tips created by
%pruning???? -HLE


%First, find all the terminal vertices (tips) and their associated edges
[iTipVert,iTipEdge] = findTips(edges,nVerts);
nTip = numel(iTipVert);


%Go through every tip and check it's length
tooShort = false(nTip,1);
for j = 1:nTip    
    %Calculate the length along this edge
    edgeLength = sum(sqrt(sum(diff(edgePaths{iTipEdge(j)},1,1) .^2,2)));
    %If it's too short, mark it for pruning
    if edgeLength < p.MinLength
        tooShort(j) = true;
        if p.ShowPlots
            if ~exist('mlFig','var')
                mlFig = plotSkel(vertices,edges,edgePaths);
                title('Minimum Length Pruning')
            end            
            plot3(edgePaths{iTipEdge(j)}(:,2),edgePaths{iTipEdge(j)}(:,1),edgePaths{iTipEdge(j)}(:,3),'--xk','LineWidth',2,'MarkerSize',10) 
        end
        
    end    
end

%Remove the edges which were too short
[vertices,edges,edgePaths,~,~] = deleteEdges(vertices,edges,edgePaths,iTipEdge(tooShort));
%This may have introduced degree=2 vertices, so once again we need to
%remove these
[vertices,edges,edgePaths,nEdges,nVerts] = removeDegree2Vert(vertices,edges,edgePaths,p);


%% --------------- Image - Edge Pruning -------------- %%
%Removes any branch tips that touch the image border.

%Find tips again as they may have changed due to the previous pruning
[iTipVert,iTipEdge] = findTips(edges,nVerts);
nTip = numel(iTipVert);

%Loop through each tip and check how close to the image edge it is
tooClose = false(nTip,1);
for j = 1:nTip    
    if any([M-vertices(iTipVert(j),1) N-vertices(iTipVert(j),2), P-vertices(iTipVert(j),3)] < p.ImageEdgeDist) || ...
        any(vertices(iTipVert(j),:) < p.ImageEdgeDist)
        
        tooClose(j) = true;
        
        if p.ShowPlots
            if ~exist('ieFig','var')
                ieFig = plotSkel(vertices,edges,edgePaths); %#ok<NASGU>
                title('Image Edge Pruning')
            end            
            plot3(edgePaths{iTipEdge(j)}(:,2),edgePaths{iTipEdge(j)}(:,1),edgePaths{iTipEdge(j)}(:,3),'--xk','LineWidth',2,'MarkerSize',10) 
        end        
    end        
end

%Remove the edges which were too close to the image edge
[vertices,edges,edgePaths,~,~] = deleteEdges(vertices,edges,edgePaths,iTipEdge(tooClose));
%This may have introduced degree=2 vertices, so once again we need to
%remove these
[vertices,edges,edgePaths,nEdges,nVerts] = removeDegree2Vert(vertices,edges,edgePaths,p);



%% ----------- Aspect Ratio Pruning, Branch/Body assignment ------------ %%
%Assigns edges to either the body or branches of the skeleton based on
%geometric criteria


%Get the depth of each edge within the mask. This is used as an
%approximation of the radius of that portion of the mask.
edDepth = edgeDepths(edgePaths,mask);
%we will also need the 


%Yet again calculate the degree of each vertex
vDegree = graphVertDegree(edges,nVert);

%Go through each edge and prune or label it
edgeLabels = nan(nEdges,1);
for j = 1:nEdges
    
    %Check the degree of the vertices of this edge
    edVertDeg = vDegree(edges(j,:));
                
    if any(edVertDeg == 1)

        
        %Get the length of each segment of this edge
        segLength = sqrt(sum(diff(edgePaths{j},1,1) .^2,2));

        %Get the index of the degree-1 vertex
%        iDeg1 = edges(j,edVertDeg==1);
%         %Find which end of the edgepath connects to this vertex
%         dStart = sqrt(sum((vertices(iDeg1,:)-edgePaths{j}(1,:)) .^2));
%         dEnd = sqrt(sum((vertices(iDeg1,:)-edgePaths{j}(end,:)) .^2));
%                 
%         %Flip the edge-Path so that it starts next to the deg=1 vertex.
%         %This simplifies the commands below and will be helpful later
%         if dEnd < dStart
%             edgePaths{j} = edgePaths{j}(end:-1:1,:);
%         end                
%         
%                 
%         if any(edDepth{j} > p.MaxRadius)
%             
%             %Use only the aspect ratio of the 
%         
%         
%             
%         else
%             %Calculate the average aspect ratio of the whole branch
%             aspRat = sum(segLength) / mean(edDepth{j});                               
%         end

        aspRat = sum(segLength) / mean(edDepth{j});
                           
        %assign to body/branches based on this ratio
        if aspRat > p.MinAspectRatio
            edgeLabels(j) = 1;%Label it as a branch
        else
            edgeLabels(j) = 2;%Label it as cell body
        end
                        
    else
        %If we are not at a branch tip, we classify only by radius
        if mean(edDepth{j}) < p.MaxRadius
            edgeLabels(j) = 1; %Label it as a branch
        else
            edgeLabels(j) = 2; %Label this as part of the cell body
        end
    end
    
    if p.ShowPlots
        if ~exist('apFig','var')
            apFig = plotSkel(vertices,edges,edgePaths); %#ok<NASGU>
            patch(maskProp.SmoothedSurface,'FaceAlpha',.1,'EdgeColor','none');
            light;
            title('Aspect ratio-based assignment. Black=Branches, Red=body')
        end
        if edgeLabels(j) == 1
            plot3(edgePaths{j}(:,2),edgePaths{j}(:,1),edgePaths{j}(:,3),'--k','LineWidth',3,'MarkerSize',10)
        else
            plot3(edgePaths{j}(:,2),edgePaths{j}(:,1),edgePaths{j}(:,3),'--r','LineWidth',3,'MarkerSize',10)
        end
        
    end    
    
end


% %Determine which vertices are end-points (branch tips). These are vertices
% %which only connect to one edge.
% nEdgesPerVert = zeros(nVert,1);
% for j = 1:nVert
%     %Count the number of edges which this vertex is connected to
%     nEdgesPerVert(j) = nnz(edges == j);       
% end
% 
% %Find endpoints (they only connect to one edge)
% iEndPt = find(nEdgesPerVert == 1);
% nEndPts = numel(iEndPt);
% 
% %Find edge which this endpoint is the end of, and the vertex that starts it
% iStartPt = zeros(nEndPts,1);
% iEdge = zeros(nEndPts,1);
% for j = 1:nEndPts
% 
%     %Get the index of the edge for this point.
%     tmp = find(arrayfun(@(x)(any(edges(x,:) == iEndPt(j))),1:nEdges) );
%     
%     if ~isempty(tmp)%Make sure it isn't a spur
%         iEdge(j) = tmp;
%         if edges(iEdge(j),1) == iEndPt(j)
%             iStartPt(j) = edges(iEdge(j),2);
%         else
%             iStartPt(j) = edges(iEdge(j),1);
%         end        
%     else
%        iStartPt(j) = iEndPt(j); 
%     end
% end
% 
% if showPlots
%     %Show the skeleton prior to pruning
%     fsFigure(.75);
%     hold on
%     patch(maskSurf,'EdgeColor','none','EdgeAlpha',.1,'FaceAlpha',.2)
%     axis vis3d,axis equal,light
%     cellfun(@(x)(plot3(x(:,2),x(:,1),x(:,3),'k','LineWidth',2)),edgePaths(cellfun(@(x)(~isempty(x)),edgePaths)));    
%     arrayfun(@(x)(text(vertices(iEndPt(x),2),vertices(iEndPt(x),1),vertices(iEndPt(x),3),num2str(x),'color','r')),1:nEndPts);
%     title('Original, un-pruned skeleton')
% end
% 
% 
% %We have curvature data for the faces, so we want their positions also.
% %Average of the vertex locations gives us the barycenter of each face. 
% facePos = zeros(nFaces,3);
% facePos(:,1) = arrayfun(@(x)(mean(maskSurf.vertices(maskSurf.faces(x,:),1),1)),1:nFaces);
% facePos(:,2) = arrayfun(@(x)(mean(maskSurf.vertices(maskSurf.faces(x,:),2),1)),1:nFaces);
% facePos(:,3) = arrayfun(@(x)(mean(maskSurf.vertices(maskSurf.faces(x,:),3),1)),1:nFaces);
% 
% 
% %Find "depth" at each point within mask via distance transform
% distX = bwdist(~maskIn);
% 
% 
% %% ------------ Pruning ---------------- %%
% 
% 
% %Get depth at beginning and end of branches.
% endDepth = arrayfun(@(x)(distX(round(vertices(x,1)),round(vertices(x,2)),round(vertices(x,3)) )),iEndPt);
% startDepth = arrayfun(@(x)(distX(round(vertices(x,1)),round(vertices(x,2)),round(vertices(x,3)) )),iStartPt);
% 
% %Calculate local surface curvature near every endpoint
% meanTipCurvature = zeros(nEndPts,1);
% gaussTipCurvature = zeros(nEndPts,1);
% k1TipCurvature = zeros(nEndPts,1);
% k2TipCurvature = zeros(nEndPts,1);
% isGoodEP = true(nEndPts,1);
% branchRadius = cell(nEndPts,1);
% 
% fitCoef = zeros(2,nEndPts);
% 
% for j = 1:nEndPts
% 
%     %Find the closest vertex to this point
%     [~,iClosest] = min(arrayfun(@(x)(sqrt(sum((maskSurf.vertices(x,[2 1 3]) - ...
%         vertices(iEndPt(j),:)).^2))),1:nSurfVert));%Vertices are in matrix coord, so we need to rearrange
% 
%     %Find which faces are near this on the surface        
%     [~,isCloseEnough] = adjacentMeshElements(maskSurf,iClosest,(curvSampD+endDepth(j)));
%     %De-cell the indices
%     isCloseEnough = isCloseEnough{1};
% 
%     %Get average mean and gaussian curvature   
%     if nnz(isCloseEnough) > 0
%         meanTipCurvature(j) = mean(H(isCloseEnough));
%         gaussTipCurvature(j) = mean(K(isCloseEnough));
%         k1TipCurvature(j) = mean(k1(isCloseEnough));
%         k2TipCurvature(j) = mean(k2(isCloseEnough));
%          
%     else
%        %If we couldn't find any surface points within the search radius, we
%        %give the branch -Inf curvature so it will pass maximum curvature
%        %criteria.
%        meanTipCurvature(j) = -Inf;
%        gaussTipCurvature(j) = -Inf;
%        k1TipCurvature(j) = -Inf;
%        k2TipCurvature(j) = -Inf;
%     end
%     
%     if showPlots            
%         if ~isempty(isCloseEnough)
%             plot3(facePos(isCloseEnough,1),facePos(isCloseEnough,2),facePos(isCloseEnough,3),'b.')
%         end
%         plot3(maskSurf.vertices(iClosest,1),maskSurf.vertices(iClosest,2),maskSurf.vertices(iClosest,3),'gx')
%     end
%     
% 
%     %Sample the distance transform along the branch corresponding to this
%     %endpoint
%     if iEdge(j) ~= 0 %First make sure it's not a spur 
%         nEdgePts = size(edgePaths{iEdge(j)},1);
%         branchRadius{j} = arrayfun(@(x)(distX(round(edgePaths{iEdge(j)}(x,1)),...
%                                               round(edgePaths{iEdge(j)}(x,2)),...
%                                               round(edgePaths{iEdge(j)}(x,3)))),...
%                                               1:nEdgePts);
% 
%         if nEdgePts >= 3
% 
%             [fitCoef(:,j),tmp] = robustfit(1:nEdgePts,branchRadius{j});
% 
%             %Calculate R^2 for this fit
%             lFun = @(x)(x * fitCoef(2,j) + fitCoef(1,j));
%             tmp.Rsquared = 1 - sum((branchRadius{j} - lFun(1:nEdgePts)) .^2) / ...
%                                sum((branchRadius{j} - mean(branchRadius{j})) .^2);
%                            
%             if j == 1
%                 %Initialize all the structures to these fields.
%                 fitStats = repmat(tmp,1,nEndPts);
%             end
%             fitStats(j) = tmp; %Allow extra field for Rsquared                                    
% 
%             if fitStats(j).p(2) < .1 && fitStats(j).Rsquared > .7 && abs(fitCoef(2,j) > .25)
%                 isGoodEP(j) = false;
% 
%             end
%         end
%     end
% end
% 
% %Find endpoints for which the surface mean curvature is convex.
% isGoodEP(meanTipCurvature > (mcThresh/4)) = false;
% 
% 
% %Check the depth of the start-point for the branches with endpoints. If
% %these are too deep, this is not a "real" branch
% %TEMP - NOT DOING THIS CHECK YET - HLE
% 
% 
% 
% %Remove branches whose start point is too deep.
% isGoodEP(startDepth > maxRad) = false;
% 
% 
% %% ---------- Output ---------- %%
% 
% 
% %Check if there are any vertices which now have no edges and delete them
% nVert = size(vertices,1);
% vDeg = graphVertDegree(edges,nVert);
% 
% 
% %Return the pruned skeleton graph
% vertices = vertices(iEndPt(isGoodEP),:);
% edgePaths(iEdge(~isGoodEP & iEdge > 0)) = [];
% edges(iEdge(~isGoodEP & iEdge > 0),:) = [];
% 
% if showPlots
%     fsFigure(.75);    
%     hold on
%     patch(maskSurf,'FaceColor','flat','EdgeColor','none','FaceVertexCData',K,'AmbientStrength',.75,'FaceAlpha',.3)
%     caxis([nanmean(K)-2*nanstd(K) nanmean(K)+2*nanstd(K)])
%     hold on                 
%     plot3(vertices(:,2),vertices(:,1),vertices(:,3),'or','MarkerSize',10);
%     cellfun(@(x)(plot3(x(:,2),x(:,1),x(:,3),'k','LineWidth',2)),edgePaths(cellfun(@(x)(~isempty(x)),edgePaths)));    
%     %arrayfun(@(x)(spy3d(vertices == x,'or','MarkerSize',15)),1:nVerts)    
%     light
%     axis image,axis vis3d    
%     view(3)
%     title('Final, Pruned Skeleton Graph with Mask Surface')           
% end

function [vertices,edges,edgePaths,nEdges,nVerts] = deleteEdges(vertices,edges,edgePaths,iEdge)



%Remove the edge and edgePath
edges(iEdge,:) = [];
edgePaths(iEdge) = [];

%Return the new edge and vertex numbers
nEdges = size(edges,1);
nVerts = size(vertices,1);

function [vertices,edges,edgePaths,nEdges,nVerts] = removeDegree2Vert(vertices,edges,edgePaths,p)

%First calcualte the degree of each vertex
vDeg = graphVertDegree(edges,size(vertices,1));

iDeg2 = find(vDeg==2);
nDeg2 = numel(iDeg2);
iPrune = nan(nDeg2,1);


for j = 1:nDeg2
    
    %Get the indices of the two edges for this vertex
    iEdges = find(any(bsxfun(@eq,edges,iDeg2(j)),2));
    
    %Merge the second edge with the first one
    edgePaths{iEdges(1)} = mergeEdgePaths(edgePaths(iEdges));        
    
    %Mark the second edge for deletion    
    iPrune(j) = iEdges(2);            
    %Get the other vertex the first edge connects to
    iVert1  = edges(iEdges(1),edges(iEdges(1),:)~=iDeg2(j));
    %... and the other vertex the merged edge connects to
    iVert2 = edges(iEdges(2),edges(iEdges(2),:)~=iDeg2(j));    
    %Update the edge specification to reflect the merge and removal of the
    %degree 2 vertex
    edges(iEdges(1),:) = [iVert1 iVert2];
    edges(iEdges(2),:) = 0;%Un-assign these vertices to prevent this edge from being merged later
    
    if p.ShowPlots
        
        if ~exist('d2Fig','var')
            d2Fig = plotSkel(vertices,edges,edgePaths); %#ok<NASGU>
            title('Degree 2 vertex removal')
        end        
        plot3(edgePaths{iEdges(1)}(:,2),edgePaths{iEdges(1)}(:,1),edgePaths{iEdges(1)}(:,3),'--k','LineWidth',2)
        plot3(vertices(iDeg2(j),2),vertices(iDeg2(j),1),vertices(iDeg2(j),3),'rx','MarkerSize',20)                
        
    end
   
    
end

%Remove the edges from the skeleton graph which were merged. This will also
%remove the vertices which are now unconnected due to the merge
[vertices,edges,edgePaths,nEdges,nVerts] = deleteEdges(vertices,edges,edgePaths,iPrune);

function edgePath = mergeEdgePaths(edgePaths)

%If the two edges to merge are further apart than this distance, a warning
%will be produced
maxMerge = 10;


iPts  = [ 1                        1;
          1                     size(edgePaths{2},1);
          size(edgePaths{1},1)     1;
          size(edgePaths{1},1)  size(edgePaths{2},1)];


%Calc distance between each start/end point pair
seDist = nan(4,1);
for j = 1:4
    seDist(j) = sqrt(sum((edgePaths{1}(iPts(j,1),:) - edgePaths{2}(iPts(j,2),:)) .^2));
end

[minD,iMinD] = min(seDist);

if minD > maxMerge
    warning('pruneSkel:mergeDistance','Two edges were merged, but the closest distance between them was %d',minD)
end

%Combine the two edges, preserving the directionality of the first edge
switch iMinD    
    
    case 1
        %The start points are closest - first one goes backwards, second
        %forwards
        edgePath = vertcat(edgePaths{1}(end:-1:1,:),edgePaths{2});
        
    case 2
        %Start of 1 is closest to end of 2 - both go backwards
        edgePath = vertcat(edgePaths{1}(end:-1:1,:),edgePaths{2}(end:-1:1,:));
        
    case 3
        %End of 1 is closest to start of 2 - both go forwards
        edgePath = vertcat(edgePaths{:});
        
    case 4 
        %End points are closest - second one goes backwards
        edgePath = vertcat(edgePaths{1},edgePaths{2}(end:-1:1,:));
        
end

function [iTipVert,iTipEdge] = findTips(edges,nVerts)

%First, get the degree of each vertex and find tips (deg = 1)
vDegree = graphVertDegree(edges,nVerts);
iTipVert = find(vDegree == 1);
nTip = numel(iTipVert);
iTipEdge = zeros(nTip,1);
%Now get the edge associated with each tip
for j = 1:nTip    
    iTipEdge(j) = find(any(bsxfun(@eq,edges,iTipVert(j)),2));    
end


    









