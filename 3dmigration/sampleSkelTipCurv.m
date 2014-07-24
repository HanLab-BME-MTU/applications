function tipCurv = sampleSkelTipCurv(vertices,edges,edgePaths,maskProp,mask,varargin)




% TEMP - UNDER CONSTRUCTION!!@!!!!!!!
%ADD HELP!!

%% ----------------------- Parameters ----------------- %%

projDist = 10;%Initial projection distance - how far from an edge to project onto the surface to sample curvature.
maxPrj = 2;   %Maximum multiple of the projection distance of a tip from the surface for it to be 
              %sampled, in the direction of the tip's edge. Tips which are
              %further than this will have NaN as their sampled curvatures.

%% ------------------------ Input ------------------------ %%

%Parse all the inputs
ip = inputParser;
ip.FunctionName = mfilename;

ip.addRequired('vertices',@(x)(ndims(x) == 2 && size(x,2)==3 && size(x,1)>1))
ip.addRequired('edges',@(x)(ndims(x) == 2 && size(x,2)==2 && size(x,1)>1))
ip.addRequired('edgePaths',@(x)(iscell(x) && numel(x) > 1));
ip.addRequired('maskProp',@(x)(isstruct(x) && numel(x) == 1));
ip.addRequired('mask',@(x)(ndims(x) == 3 && islogical(x)));
ip.addParamValue('ShowPlots',false,(@(x)(numel(x)==1)));
ip.addParamValue('CurvSampRad',5,(@(x)(numel(x) == 1 && x >= 1))); 

ip.parse(vertices,edges,edgePaths,maskProp,mask,varargin{:});
p = ip.Results;


%% ----------------------- Init -------------------------- %%

nVert = size(vertices,1);

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
    

    %Find the mask surface vertex associated with this tip
    
    %First, find the edge connecting to this tip.        
    isCurrEdge = any(bsxfun(@eq,edges,iTips(j)),2);
    
    %Find which end of the edgepath connects to this vertex
    dStart = sqrt(sum((vertices(iTips(j),:)-edgePaths{isCurrEdge}(1,:)) .^2));
    dEnd = sqrt(sum((vertices(iTips(j),:)-edgePaths{isCurrEdge}(end,:)) .^2));

    %Flip the edge-Path so that it ends next to the tip
    if dStart < dEnd
        edgePaths{isCurrEdge} = edgePaths{isCurrEdge}(end:-1:1,:);
    end                
    
    
    %Project this edge onto the mask surface to find the surface point
    %associated with this tip. We keep projecting the edge until we hit the surface.
    %TEMP - you can speed this up - don't re-interp the old points! (but it
    %happens only rarely so...)
    surfPoint = [];
    nProj = 1;
    while isempty(surfPoint) && nProj <= maxPrj
    
        %Get the average direction of the last portion of the edge path        
        if size(edgePaths{isCurrEdge},1) > 1
            [~,edgeProj] = curveDirection3D(edgePaths{isCurrEdge},4,projDist*nProj);    
        else
            %If the edge is too short, just use the line connecting the two
            %vertices
            [~,edgeProj] = curveDirection3D(vertices(edges(isCurrEdge,:),:),4,projDist*nProj);    
            
        end
        %Use this to sample to mask to find the point where it hits the mask
        %surface
        maskSamp = interp3(single(mask),edgeProj(:,2),edgeProj(:,1),edgeProj(:,3));
        maskSamp(isnan(maskSamp)) = 0;%If we hit the image border, stop there.
        surfPoint = edgeProj(find(maskSamp == 0,1,'first'),:);
        nProj = nProj + 1;
        
    end
    
    if ~isempty(surfPoint)
    
        %Find which faces are near this on the surface        
        [~,isCloseEnough] = adjacentMeshElements(maskProp.SmoothedSurface,surfPoint([2 1 3]),p.CurvSampRad);
        %De-cell the indices
        if ~isempty(isCloseEnough)
            isCloseEnough = isCloseEnough{1};       
        end
    else
        isCloseEnough = 0;
    end
            
    %Get average curvatures for this tip    
    if nnz(isCloseEnough) > 0
        
        tipCurv.meanTipCurvature(iTips(j)) = robustMean(maskProp.MeanCurvature(isCloseEnough));
        tipCurv.gaussTipCurvature(iTips(j)) = robustMean(maskProp.GaussianCurvature(isCloseEnough));
        tipCurv.k1TipCurvature(iTips(j)) = robustMean(maskProp.CurvaturePC1(isCloseEnough));
        tipCurv.k2TipCurvature(iTips(j)) = robustMean(maskProp.CurvaturePC2(isCloseEnough));
        tipCurv.diffk1k2TipCurvature(iTips(j)) = robustMean(real(maskProp.CurvaturePC1(isCloseEnough)) - ...
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
        if nnz(isCloseEnough)>0
            %Plot one vertex of every sampled face
            plot3(maskProp.SmoothedSurface.vertices(maskProp.SmoothedSurface.faces(isCloseEnough,1),1),...
                  maskProp.SmoothedSurface.vertices(maskProp.SmoothedSurface.faces(isCloseEnough,1),2),...
                  maskProp.SmoothedSurface.vertices(maskProp.SmoothedSurface.faces(isCloseEnough,1),3),'b.')
        end
        if ~isempty(surfPoint)
            plot3(surfPoint(2),surfPoint(1),surfPoint(3),'gx')
        end
    end
        
end


