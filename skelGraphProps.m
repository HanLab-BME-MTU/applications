function grphStats= skelGraphProps(vertices,edges,edgePaths,mask,maskProp,varargin)
%SKELGRAPHSTATS calculates various properties & statistics of the input skeleton graph
%
% grphStats = skelGraphStats(vertices,edges,edgePaths,mask,maskProp)
% grphStats = skelGraphStats(...,'OptionName1',optionValue1,...)
%
% This function calculates various geometric and topological properties for
% each element of the input skeleton graph using the input mask and mask
% geometry analysis. For a list of calculated properties, see "Output"
%
% 
% Input:
%
%   vertices,edges,edgePaths - The graph form of the skeleton, in the
%   format used by skel2graph.m
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
%       ('CurvSampRadius'->scalar >= 1) The radius over which to
%       sample the local surface curvature of a branch tip. The tip
%       curvature will be equal to the average curvature of all surface
%       points within this distance of the skeleton branch tip. Optional.
%       Default is 2.
%
%       ('ShowPlots'->true/false) If true, figures will be displayed
%       showing various calculated properties
%   
%
% Output:
%
%   grphStats - A structure containing the calculated skeleton properties.
%       Fields include:
%           -vertexDegree - The number of edges connecting to each vertex           
%           ****TEMP - EXPAND!!!!! ******
%
% Hunter Elliott
% 5/2011
%

%% ------------------------ Input ------------------------ %%

%Parse all the inputs
ip = inputParser;
ip.FunctionName = mfilename;

ip.addRequired('vertices',@(x)(ndims(x) == 2 && size(x,2)==3 && size(x,1)>1))
ip.addRequired('edges',@(x)(ndims(x) == 2 && size(x,2)==2 && size(x,1)>1))
ip.addRequired('edgePaths',@(x)(iscell(x) && numel(x) > 1));
ip.addRequired('mask',@(x)(ndims(x) == 3 && islogical(x)));
ip.addRequired('maskProp',@(x)(isstruct(x) && numel(x) == 1));
ip.addParamValue('ShowPlots',false,(@(x)(numel(x)==1)));
ip.addParamValue('CurvSampRad',2,(@(x)(numel(x) == 1 && x >= 1))); 

ip.parse(vertices,edges,edgePaths,mask,maskProp,varargin{:});
p = ip.Results;


%% ----------------------- Init -------------------------- %%

nVert = size(vertices,1);
nEdges = size(edges,1);
grphStats.vertexDegree = zeros(max(edges(:)),1);
nSurfVert = size(maskProp.SmoothedSurface.vertices,1);

%% -------------------- Property Calculation ------------ %%

% ---- Vertex Degree ----- %

for j = 1:nEdges    
    %Make sure it's not a spur
    if ~any(edges(j,:)==0)
        %Add this edge to the degree count of each vertex it connects
        grphStats.vertexDegree(edges(j,:)) = grphStats.vertexDegree(edges(j,:)) + 1;
    end
end

if p.ShowPlots    
    fsFigure(.75);
    hold on,axis vis3d,axis equal
    cellfun(@(x)(plot3(x(:,1),x(:,2),x(:,3),'k','LineWidth',2)),edgePaths(cellfun(@(x)(~isempty(x)),edgePaths)));            
    for j = 1:nVert        
        plot3(vertices(j,1),vertices(j,2),vertices(j,3),'or','MarkerSize',15)
        text(vertices(j,1),vertices(j,2),vertices(j,3),num2str(grphStats.vertexDegree(j)),'Color','b')
    end
    title('Skeleton Graph Vertex Degree')
end

% ---- Edge Multiplicity ---- %

%Calculate how many edges connect each pair of vertices
grphStats.edgeMultiplicity = nan(nEdges,1);
for j = 1:nEdges            
    grphStats.edgeMultiplicity(j) = nnz(all(bsxfun(@eq,edges,edges(j,:)),2));    
end



% ---- Tip Curvature ---- %

%Find "tips" - vertices with degree = 1
iTips = find(grphStats.vertexDegree == 1);
nTips = numel(iTips);

grphStats.meanTipCurvature = nan(nVert,1);
grphStats.gaussTipCurvature = nan(nVert,1);
grphStats.k1TipCurvature = nan(nVert,1);
grphStats.k2TipCurvature = nan(nVert,1);
grphStats.diffk1k2TipCurvature = nan(nVert,1);

for j = 1:nTips
    
    %Find the mask surface vertex closest to this tip
    %TEMP - if far from surface, project branch direction instead of using
    %closest point??? - HLE
    [~,iClosest] = min(arrayfun(@(x)(sqrt(sum((maskProp.SmoothedSurface.vertices(x,[2 1 3]) - ...
        vertices(iTips(j),:)).^2))),1:nSurfVert));%Vertices are in matrix coord, so we need to rearrange
    
    %Find which faces are near this on the surface        
    [~,isCloseEnough] = adjacentMeshElements(maskProp.SmoothedSurface,iClosest,p.CurvSampRad);
    %De-cell the indices
    isCloseEnough = isCloseEnough{1};
    
    %Get average curvatures for this tip    
    if nnz(isCloseEnough) > 0
        
        grphStats.meanTipCurvature(iTips(j)) = nanmean(maskProp.MeanCurvature(isCloseEnough));
        grphStats.gaussTipCurvature(iTips(j)) = nanmean(maskProp.GaussianCurvature(isCloseEnough));
        grphStats.k1TipCurvature(iTips(j)) = nanmean(maskProp.CurvaturePC1(isCloseEnough));
        grphStats.k2TipCurvature(iTips(j)) = nanmean(maskProp.CurvaturePC2(isCloseEnough));
        grphStats.diffk1k2TipCurvature(iTips(j)) = nanmean(real(maskProp.CurvaturePC1(isCloseEnough)) - ...
                                                           real(maskProp.CurvaturePC2(isCloseEnough)));                         
    end
    
    if p.ShowPlots            
        if j == 1
            fsFigure(.75);
            hold on
            patch(maskProp.SmoothedSurface,'FaceAlpha',.3,'EdgeColor','none');
            axis vis3d, axis equal,view(3),light
            cellfun(@(x)(plot3(x(:,2),x(:,1),x(:,3),'k','LineWidth',2)),edgePaths(cellfun(@(x)(~isempty(x)),edgePaths)));                        
            arrayfun(@(x)(text(vertices(x,2),vertices(x,1),vertices(x,3),num2str(x),'Color','r')),1:nVert);
        end
        if ~isempty(isCloseEnough)
            %Plot one vertex of every sampled face
            plot3(maskProp.SmoothedSurface.vertices(maskProp.SmoothedSurface.faces(isCloseEnough,1),1),...
                  maskProp.SmoothedSurface.vertices(maskProp.SmoothedSurface.faces(isCloseEnough,1),2),...
                  maskProp.SmoothedSurface.vertices(maskProp.SmoothedSurface.faces(isCloseEnough,1),3),'b.')
        end
        plot3(maskProp.SmoothedSurface.vertices(iClosest,1),maskProp.SmoothedSurface.vertices(iClosest,2),maskProp.SmoothedSurface.vertices(iClosest,3),'gx')
    end
        
end

% ---- Cumulative Length & Per-Point Depth  ---- %



%Get the distance transform for calculating depth of branch points
distX = bwdist(~mask);

grphStats.cumEdgeLength = cell(nEdges,1);
grphStats.edgeDepth = cell(nEdges,1);
for j = 1:nEdges
    
    %Calculate the cumulative length along this edge
    grphStats.cumEdgeLength{j} = cumsum(sqrt(sum(diff(edgePaths{j},1,1) .^2,2)));        
    
    
    %Get the depth within the mask of each point on this edge
    nEdgePts = size(edgePaths{j},1);
    grphStats.edgeDepth{j} = arrayfun(@(x)(distX(round(edgePaths{j}(x,1)),...
                                                 round(edgePaths{j}(x,2)),...
                                                 round(edgePaths{j}(x,3)))),...
                                                 1:nEdgePts)';

    if p.ShowPlots
        if j == 1
            fsFigure(.75);
            hold on
            patch(maskProp.SmoothedSurface,'FaceAlpha',.2,'EdgeColor','none');
            axis vis3d, axis equal,view(3),light
            cols = jet(nEdges);
        end
        text(mean(edgePaths{j}(:,2)),mean(edgePaths{j}(:,1)),mean(edgePaths{j}(:,3)),num2str(j),'color',cols(j,:));
        plot3(edgePaths{j}(:,2),edgePaths{j}(:,1),edgePaths{j}(:,3),'color',cols(j,:),'LineWidth',2);
    end                           
    
end



