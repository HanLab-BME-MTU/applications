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
curvSlack = 3;%Tip curvature must be below maxRad by at least this factor 
              %for a tip to be pruned. This is because a branch with radius
              %<= MaxRadius can still have a fairly flat tip depending on
              %it's geometry, so we give some leeway

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
ip.addParamValue('MinAspectRatio',4,(@(x)(numel(x) == 1 && x >= 1))); 
ip.addParamValue('MinLength',5,(@(x)(numel(x) == 1 && x >= 1))); 
ip.addParamValue('ImageEdgeDist',3,(@(x)(numel(x) == 1 && x >= 1)));
ip.addParamValue('CurvSampRad',7,(@(x)(numel(x) == 1 && x >= 1))); 

ip.parse(vertices,edges,edgePaths,mask,maskProp,varargin{:});
p = ip.Results;


%TEMP - MORE CHECKING OF REQUIRED INPUTS???!!-HLE


%% -------------------------- Init ------------------------------%%

%Get radius-of-curvature thresholds from the maximum radius parameter
gcThresh = 1/p.MaxRadius^2;%Maximum gaissian curvature value. This corresponds to
                      %the gaussian curvature of a sphere of radius maxRad
                         
mcThresh = -1/p.MaxRadius; %Maximum mean curvature value. This corresponds to
                      %the mean curvature of a sphere of radius maxRad

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

                plot3(edgePaths{iSharing(2)}(:,2),edgePaths{iSharing(2)}(:,1),edgePaths{iSharing(2)}(:,3),'rx','MarkerSize',20)
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

%% ---------------- Surface Curvature Pruning ------------------ %%
%Removes branches whose tips touch the surface at areas of very low
%curvature


%Sample the curvature of the mask surface adjacent to each tip
tipCurv = sampleSkelTipCurv(vertices,edges,edgePaths,maskProp,'ShowPlots',p.ShowPlots,'CurvSampRad',p.CurvSampRad);
%Get the indices of tips again, as this may have changed due to last
%pruning
[iTipVert,iTipEdge] = findTips(edges,nVerts);
nTip = numel(iTipVert);
tooFlat = false(nTip,1);

%Go through each tip and check the surface curvature
for j = 1:nTip
    
    %Check both mean and gaussian curvatures, since the curvatures are
    %approximate and slightly noisy.
    if tipCurv.meanTipCurvature(iTipVert(j)) > (mcThresh/curvSlack) && ... %use > for mean, because or system gives negative values for convex areas
        tipCurv.gaussTipCurvature(iTipVert(j)) < (gcThresh/curvSlack)
        tooFlat(j) = true;
        if p.ShowPlots
            plot3(edgePaths{iTipEdge(j)}(:,2),edgePaths{iTipEdge(j)}(:,1),edgePaths{iTipEdge(j)}(:,3),'--xr','LineWidth',2,'MarkerSize',10) 
        end
    elseif p.ShowPlots
        plot3(edgePaths{iTipEdge(j)}(:,2),edgePaths{iTipEdge(j)}(:,1),edgePaths{iTipEdge(j)}(:,3),'--g','LineWidth',2,'MarkerSize',10) 
    end
end

%Remove the edges marked for pruning
%Remove the edges which were too close to the image edge
[vertices,edges,edgePaths,~,~] = deleteEdges(vertices,edges,edgePaths,iTipEdge(tooFlat));
%Once again, remove deg=2 vertices
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
        iDeg1 = find(edges(j,edVertDeg==1),1);
        %Find which end of the edgepath connects to this vertex
        dStart = sqrt(sum((vertices(iDeg1,:)-edgePaths{j}(1,:)) .^2));
        dEnd = sqrt(sum((vertices(iDeg1,:)-edgePaths{j}(end,:)) .^2));
                
        %Flip the edge-Path so that it starts next to the deg=1 vertex.
        %This simplifies the commands below and will be helpful later
        if dEnd < dStart
            edgePaths{j} = edgePaths{j}(end:-1:1,:);
        end                
        
                
        if any(edDepth{j} > p.MaxRadius)
            
            %Use only the length of the portion with radius less than
            %maxradius. This is sort of a cheap trick, but it seems to
            %work...
            iLastGoodRad = find(edDepth{j}<p.MaxRadius,1,'Last');
            aspRat = sum(segLength(1:max(iLastGoodRad-1,1))) / mean(edDepth{j});
            
        else
            %Calculate the average aspect ratio of the whole branch
            aspRat = sum(segLength) / mean(edDepth{j});                               
        end
        
                           
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


    









