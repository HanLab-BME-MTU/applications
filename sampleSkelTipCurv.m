function tipCurv = sampleSkelTipCurv(vertices,edges,edgePaths,maskProp,varargin)




% TEMP - UNDER CONSTRUCTION!!@!!!!!!!
%ADD HELP!!


%% ------------------------ Input ------------------------ %%

%Parse all the inputs
ip = inputParser;
ip.FunctionName = mfilename;

ip.addRequired('vertices',@(x)(ndims(x) == 2 && size(x,2)==3 && size(x,1)>1))
ip.addRequired('edges',@(x)(ndims(x) == 2 && size(x,2)==2 && size(x,1)>1))
ip.addRequired('edgePaths',@(x)(iscell(x) && numel(x) > 1));
%ip.addRequired('mask',@(x)(ndims(x) == 3 && islogical(x)));
ip.addRequired('maskProp',@(x)(isstruct(x) && numel(x) == 1));
ip.addParamValue('ShowPlots',false,(@(x)(numel(x)==1)));
ip.addParamValue('CurvSampRad',5,(@(x)(numel(x) == 1 && x >= 1))); 

ip.parse(vertices,edges,edgePaths,maskProp,varargin{:});
p = ip.Results;


%% ----------------------- Init -------------------------- %%

nVert = size(vertices,1);
nSurfVert = size(maskProp.SmoothedSurface.vertices,1);

%Get the degree of each vertex so we can find tips (vertices with deg = 1)
vDegree = graphVertDegree(edges,nVert);

iTips = find(vDegree == 1);
nTips = numel(iTips);

tipCurv.meanTipCurvature = nan(nVert,1);
tipCurv.gaussTipCurvature = nan(nVert,1);
tipCurv.k1TipCurvature = nan(nVert,1);
tipCurv.k2TipCurvature = nan(nVert,1);
tipCurv.diffk1k2TipCurvature = nan(nVert,1);


%% ------------------ Curvature Sampling --------------- %%

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
        
        tipCurv.meanTipCurvature(iTips(j)) = nanmean(maskProp.MeanCurvature(isCloseEnough));
        tipCurv.gaussTipCurvature(iTips(j)) = nanmean(maskProp.GaussianCurvature(isCloseEnough));
        tipCurv.k1TipCurvature(iTips(j)) = nanmean(maskProp.CurvaturePC1(isCloseEnough));
        tipCurv.k2TipCurvature(iTips(j)) = nanmean(maskProp.CurvaturePC2(isCloseEnough));
        tipCurv.diffk1k2TipCurvature(iTips(j)) = nanmean(real(maskProp.CurvaturePC1(isCloseEnough)) - ...
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




