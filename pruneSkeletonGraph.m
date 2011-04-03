function [vertices,edges,edgePaths] = pruneSkeletonGraph(vertices,edges,edgePaths,maskIn,maskProps,varargin)
%PRUNESKELETONGRAPH prunes the input skeleton graph by geometric criteria
%
% [vertices,edges,edgePaths] = pruneSkeletonGraph(vertices,edges,edgePaths,maskIn,maskProps);
% [vertices,edges,edgePaths] = pruneSkeletonGraph(vertices,edges,edgePaths,maskIn,maskProps,'OptionName1',optionValue1,'OptionName2',optionValue2,...);
%
%   This function removes branches from the input skeleton based on
%   tests of several geometric criteria. These criteria may be altered or
%   removed by the user. For details of the criteria, see the option
%   descriptions below. There must be only one objec in the skeleton, the
%   mask and the maskProps structure.
%
% Input:
%
%   vertices,edges,edgePaths - The graph form of the skeleton, as returned
%   by skel2graph.m
%
%   maskIn - The original mask from which the skeleton was derived.
%
%   maskProps - The geometric properties of the mask, as returned by
%   analyze3DMaskGeometry.m
%
%   'OptionName',optionValue - A string with an option name followed by the
%   value for that option.
% 
%   Possible Option Names:
%
%       ('OptionName' -> possible values)
%
%       ('MaxRadius'->Positive scalar) The maximum radius (in voxels) of an object to consider a
%       branch. Skeleton branches resulting from areas of the mask with
%       radius larger than this will be pruned.
%       Optional. Default is 1/10th of smallest mask dimension.
%
%       ('CurvatureSamplingRadius'->Positive scalar) The radius over which
%       to sample the local curvature of a branch tip. The average
%       curvature of all surface points within this distance of the
%       skeleton branch tip will be used to determine the radius of
%       curvature of the tip, which will then be compared to the MaxRadius.       
%       Optional. Default is 2.
%
%       ('ShowPlots'->true/false) If true, figures will be displayed
%       showing the branches before and after pruning.
%   
% Output:
%
%    vertices,edges,edgePaths - The pruned skeleton graph. Same form as
%    input skeleton graph.
%
% Hunter Elliott
% 3/2011
%

%% ------------------------- Input -------------------------------%%


if nargin < 5 || isempty(vertices) || isempty(edges) || isempty(edgePaths)...
        || isempty(maskIn) || isempty(maskProps)
    error('You must input the graph skeleton (vertices,edges,edgePaths), the mask, and the mask properties!');
end


if ~islogical(maskIn) || ndims(maskIn) ~= 3
    error('The input mask must be a 3D, logical matrix!');
end

[M,N,P] = size(maskIn);

%TEMP - MORE CHECKING OF REQUIRED INPUTS???!!-HLE

[maxRad,curvSampD,showPlots] = parseInput(varargin);

if isempty(showPlots)
    showPlots = false;
end

if isempty(maxRad)
    maxRad = min([M N P])/10;
end

if isempty(curvSampD)
    curvSampD = 2;
end


%% -------------------------- Init ------------------------------%%

%Get radius-of-curvature thresholds from the maximum radius parameter

gcThresh = 1/maxRad^2;%Maximum gaissian curvature value. This corresponds to
                      %the gaussian curvature of a sphere of radius maxRad
                         
mcThresh = -1/maxRad; %Maximum mean curvature value. This corresponds to
                      %the mean curvature of a sphere of radius maxRad


nObj = numel(maskProps);%Number of objects in mask                      
                      
if nObj > 1
    error('The mask must have only one object in it!')
end

%Extract maskprop fields for readability
maskSurf = maskProps.SmoothedSurface;
H = maskProps.MeanCurvature;
K = maskProps.GaussianCurvature;
k1 = maskProps.CurvaturePC1;
k2 = maskProps.CurvaturePC2;


%Number of skeleton vertices and edges.
nEdges = size(edges,1);
nVert = size(vertices,1);

%Number of smoothed surface faces
nFaces = size(maskSurf.faces,1);



%% ------------------------ Pre-Processing ----------------------------- %%



%Determine which vertices are end-points (branch tips). These are vertices
%which only connect to one edge.
nEdgesPerVert = zeros(nVert,1);
for j = 1:nVert
    %Count the number of edges which this vertex is connected to
    nEdgesPerVert(j) = nnz(edges == j);       
end

%Find endpoints (they only connect to one edge)
iEndPt = find(nEdgesPerVert == 1 | nEdgesPerVert == 0); %The endpoints can also have 0 edges because "spurs" can exist in the skeleton - that is, branches of length = 1, where the tip IS the branch
nEndPts = numel(iEndPt);
%Find edge which this endpoint is the end of, and the vertex that starts it
iStartPt = zeros(nEndPts,1);
iEdge = zeros(nEndPts,1);
for j = 1:nEndPts

    %Get the index of the edge for this point.
    tmp = find(arrayfun(@(x)(any(edges(x,:) == iEndPt(j))),1:nEdges) );

    if ~isempty(tmp)
        iEdge(j) = tmp;
        if edges(iEdge(j),1) == iEndPt(j)
            iStartPt(j) = edges(iEdge(j),2);
        else
            iStartPt(j) = edges(iEdge(j),1);
        end        
    else
       iStartPt(j) = iEndPt(j); 
    end
end

if showPlots
    %Show the skeleton and prior to pruning
    fsFigure(.75);
    hold on
    cellfun(@(x)(plot3(x(:,2),x(:,1),x(:,3),'k')),edgePaths);    
    arrayfun(@(x)(text(vertices(iEndPt(x),2),vertices(iEndPt(x),1),vertices(iEndPt(x),3),num2str(x),'color','r')),1:nEndPts);
    title('Original, un-pruned skeleton')
end


%We have curvature data for the faces, so we want their positions also.
%Average of the vertex locations gives us the barycenter of each face. 
facePos = zeros(nFaces,3);
facePos(:,1) = arrayfun(@(x)(mean(maskSurf.vertices(maskSurf.faces(x,:),1),1)),1:nFaces);
facePos(:,2) = arrayfun(@(x)(mean(maskSurf.vertices(maskSurf.faces(x,:),2),1)),1:nFaces);
facePos(:,3) = arrayfun(@(x)(mean(maskSurf.vertices(maskSurf.faces(x,:),3),1)),1:nFaces);


%Find "depth" at each point within mask via distance transform
distX = bwdist(~maskIn);


%% --


%Get depth at beginning and end of branches.
endDepth = arrayfun(@(x)(distX(round(vertices(x,1)),round(vertices(x,2)),round(vertices(x,3)) )),iEndPt);
startDepth = arrayfun(@(x)(distX(round(vertices(x,1)),round(vertices(x,2)),round(vertices(x,3)) )),iStartPt);

%Calculate local surface curvature near every endpoint
meanTipCurvature = zeros(nEndPts,1);
gaussTipCurvature = zeros(nEndPts,1);
k1TipCurvature = zeros(nEndPts,1);
k2TipCurvature = zeros(nEndPts,1);
isGoodEP = true(nEndPts,1);
branchRadius = cell(nEndPts,1);

fitCoef = zeros(2,nEndPts);
% fitStats = struct('R',cell(nEndPts,1),...
%                   'df',cell(nEndPts,1),...
%                'normr',cell(nEndPts,1),...
%                'Rsquared',cell(nEndPts,1));

for j = 1:nEndPts

    %Calc the distance from this endpoint to every face of the surface
    %polygon
    dToSurf = arrayfun(@(x)(sqrt(sum((facePos(x,[2 1 3]) - ...
        vertices(iEndPt(j),:)).^2))),1:nFaces);%Vertices are in matrix coord, so we need to rearrange

    %Find points near this tip on the surface
    isCloseEnough = dToSurf < (curvSampD+endDepth(j));

    %Get average mean and gaussian curvature   
    if nnz(isCloseEnough) > 0
        meanTipCurvature(j) = mean(H(isCloseEnough));
        gaussTipCurvature(j) = mean(K(isCloseEnough));
        k1TipCurvature(j) = mean(k1(isCloseEnough));
        k2TipCurvature(j) = mean(k2(isCloseEnough));

         if showPlots            
             plot3(facePos(isCloseEnough,1),facePos(isCloseEnough,2),facePos(isCloseEnough,3),'b.')
         end
    else
       %If we couldn't find any surface points within the search radius, we
       %give the branch -Inf curvature so it will pass maximum curvature
       %criteria. (This is most likely because the branch was small enough
       %to be removed by the smoothing). 
       meanTipCurvature(j) = -Inf;
       gaussTipCurvature(j) = -Inf;
       k1TipCurvature(j) = -Inf;
       k2TipCurvature(j) = -Inf;
    end

    %Sample the distance transform along the branch corresponding to this
    %endpoint
    if iEdge(j) ~= 0 %First make sure it's not a spur 
        nEdgePts = size(edgePaths{iEdge(j)},1);
        branchRadius{j} = arrayfun(@(x)(distX(round(edgePaths{iEdge(j)}(x,1)),...
                                              round(edgePaths{iEdge(j)}(x,2)),...
                                              round(edgePaths{iEdge(j)}(x,3)))),...
                                              1:nEdgePts);

        if nEdgePts >= 3
%             %Fit a line to the branch series
%             [fitCoef(j,:),tmp] = polyfit(1:nEdgePts,branchRadius{j},1);                                                       
%             

            [fitCoef(:,j),tmp] = robustfit(1:nEdgePts,branchRadius{j});

            %Calculate R^2 for this fit
            lFun = @(x)(x * fitCoef(2,j) + fitCoef(1,j));
            tmp.Rsquared = 1 - sum((branchRadius{j} - lFun(1:nEdgePts)) .^2) / ...
                               sum((branchRadius{j} - mean(branchRadius{j})) .^2);
                           
            if j == 1
                %Initialize all the structures to these fields.
                fitStats = repmat(tmp,1,nEndPts);
            end
            fitStats(j) = tmp; %Allow extra field for Rsquared                                    

            if fitStats(j).p(2) < .1 && fitStats(j).Rsquared > .7 && abs(fitCoef(2,j) > .25)
                isGoodEP(j) = false;

            end
        end
    end
end

%Find endpoints for which the surface mean curvature is convex
isGoodEP(meanTipCurvature > (mcThresh/4)) = false;


%Check the depth of the start-point for the branches with endpoints. If
%these are too deep, this is not a "real" branch
%TEMP - NOT DOING THIS CHECK YET - HLE



%Remove branches whose start point is too deep.
isGoodEP(startDepth > maxRad) = false;


%% ---------- Output ---------- %%


%Return the pruned skeleton graph
vertices = vertices(iEndPt(isGoodEP),:);
edgePaths = edgePaths(iEdge(isGoodEP)>0);
edges = edges(iEdge(isGoodEP)>0);

if showPlots
    fsFigure(.75);    
    hold on
    patch(maskSurf,'FaceColor','flat','EdgeColor','none','FaceVertexCData',K,'AmbientStrength',.75,'FaceAlpha',.3)
    caxis([mean(K)-2*std(K) mean(K)+2*std(K)])
    hold on                 
    plot3(vertices(:,2),vertices(:,1),vertices(:,3),'or','MarkerSize',10);
    %arrayfun(@(x)(spy3d(vertices == x,'or','MarkerSize',15)),1:nVerts)    
    light
    axis image,axis vis3d    
    view(3)
    title('Final, Pruned Skeleton Graph with Mask Surface')           
end


function [maxRad,curvSampD,showPlots] = parseInput(argArray)

maxRad = [];
curvSampD = [];
showPlots = [];

if isempty(argArray)
    return
end

nArg = length(argArray);

%Make sure there is an even number of arguments corresponding to
%optionName/value pairs
if mod(nArg,2) ~= 0
    error('Inputs must be as optionName / value pairs!')
end

for i = 1:2:nArg
    
    
    switch argArray{i}
        
        
        case 'MaxRadius'
            
            maxRad = argArray{i+1};
            
        case 'CurvatureSamplingRadius'
            
            curvSampD = argArray{i+1};

        case 'ShowPlots'
            showPlots = argArray{i+1};                        

        otherwise

            error(['"' argArray{i} '" is not a valid option name! Please check input!'])
    end
end                       





