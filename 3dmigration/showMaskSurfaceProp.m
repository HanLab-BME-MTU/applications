function pHan = showMaskSurfaceProp(maskProp,dispType)
%SHOWMASKSURFACEPROP displays the output of analyze3DMaskGeometry
%
% showMaskSurfaceProp(maskProp)
% showMaskSurfaceProp(maskProp,displayType)
%
% Input:
%       
%   maskProp - Mask geometric property structure, as returned by
%   analyze3DMaskGeometry.m
%
%   displayType - String specifying the type of info to display. Possible
%   options are"
%
%       'gauss' - Color-codes each surface face by it's gaussian curvature.
%       This is the default
%       'mean' - Color codes by mean curvature
%       'pc1' - color codes by largest principle component of curvature
%       'pc2' - color codes by smalles principle component of curvature
%       'curv' - color codes by maximum absolute value ofcurvature component
%       'wire' - Simple wireframe surface only
%       'surf' - Simple surface only
%
%   Shows the 3D smoothed surface with the Gaussian curvature overlain on
%   it
%
% Hunter Elliott 
% 2/2012
%

if nargin < 2 || isempty(dispType)
    dispType = 'gauss';
end

switch dispType
    
    
    case 'gauss'

        pHan = patch(maskProp.SmoothedSurface,'FaceColor','flat','EdgeColor','none','FaceVertexCData',maskProp.GaussianCurvature);
        caxis([prctile(maskProp.GaussianCurvature,5) prctile(maskProp.GaussianCurvature,95)]) 

    case 'mean'

        pHan = patch(maskProp.SmoothedSurface,'FaceColor','flat','EdgeColor','none','FaceVertexCData',maskProp.MeanCurvature);
        caxis([prctile(maskProp.MeanCurvature,5) prctile(maskProp.MeanCurvature,95)])    
        
    case 'pc1'

        pHan = patch(maskProp.SmoothedSurface,'FaceColor','flat','EdgeColor','none','FaceVertexCData',maskProp.CurvaturePC1);
        caxis([prctile(maskProp.CurvaturePC1,5) prctile(maskProp.CurvaturePC1,95)]) 
        
     case 'pc2'

        pHan = patch(maskProp.SmoothedSurface,'FaceColor','flat','EdgeColor','none','FaceVertexCData',maskProp.CurvaturePC2);
        caxis([prctile(maskProp.CurvaturePC2,95) prctile(maskProp.CurvaturePC2,5)]) 
        
     case 'curv'
        
        maxCurv = max(abs(real(maskProp.CurvaturePC1)),abs(real(maskProp.CurvaturePC2)));
        pHan = patch(maskProp.SmoothedSurface,'FaceColor','flat','EdgeColor','none','FaceVertexCData',maxCurv);
        caxis([prctile(maxCurv,5) prctile(maxCurv,95)]) 
        
    case 'wire'
        
        pHan = patch(maskProp.SmoothedSurface,'FaceColor','none','EdgeColor','k');
        
    case 'surf'
        
        pHan = patch(maskProp.SmoothedSurface,'FaceColor','k','EdgeColor','none','FaceAlpha',.2);
        
    case 'LAcurv'
        
        pHan = patch(maskProp.SmoothedSurface,'FaceColor','flat','EdgeColor','none','FaceVertexCData',maskProp.locAvgCurv.LocMeanMaxAbsCurvature);
        caxis([prctile(maskProp.locAvgCurv.LocMeanMaxAbsCurvature,5) prctile(maskProp.locAvgCurv.LocMeanMaxAbsCurvature,95)]) 
        
    case 'LAgauss'
        
        pHan = patch(maskProp.SmoothedSurface,'FaceColor','flat','EdgeColor','none','FaceVertexCData',maskProp.locAvgCurv.LocMeanGaussianCurvature);
        caxis([prctile(maskProp.locAvgCurv.LocMeanGaussianCurvature,5) prctile(maskProp.locAvgCurv.LocMeanGaussianCurvature,95)]) 
        
     case 'LAmean'
        
        pHan = patch(maskProp.SmoothedSurface,'FaceColor','flat','EdgeColor','none','FaceVertexCData',maskProp.locAvgCurv.LocMeanMeanCurvature);
        caxis([prctile(maskProp.locAvgCurv.LocMeanMeanCurvature,5) prctile(maskProp.locAvgCurv.LocMeanMeanCurvature,95)]) 
        
    otherwise
        error(['"' dispType '" is not a supported display type!'])
end

set(pHan,'VertexNormals',maskProp.SurfaceNorms(:,[1 2 3]));

axis equal
light
lighting phong
view(3)
colorbar








