function showMaskSurfaceProp(maskProp,dispType)
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
%       'wire' - Simple wireframe surface only
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

        patch(maskProp.SmoothedSurface,'FaceColor','flat','EdgeColor','none','FaceVertexCData',maskProp.GaussianCurvature)
        caxis([-.4 .4])
        
    case 'wire'
        
        patch(maskProp.SmoothedSurface,'FaceColor','none','EdgeColor','k')
        
    otherwise
        error(['"' dispType '" is not a supported display type!'])
end


axis equal
light
view(3)








