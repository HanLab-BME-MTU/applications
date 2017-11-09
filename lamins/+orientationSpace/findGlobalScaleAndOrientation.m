function [ scale, orientation, scaleV, orientationV, orientationAtMaxScale ] = findGlobalScaleAndOrientation( I, varargin)
%FINDGLOBALSCALEANDORIENTATION Finds the scales and orientation using OrientationSpace
%filtering. Scans scale using Chebyshev polynomials
%
% INPUT
% I - image as a double matrix
%
% PARAMETERS (inputParser)
% AngularOrder - order to determine orientation resolution 
%                default: 5
%                See also OrientationSpaceFilter
% Scales - list of scales as a vector or length at least 2. Adjacent pairs
%          of scales will be used as intervals over which to analyze scale
%          default: 1:4
% Mask - initial mask identifying the area of interest
% TileSize - number of pixels to process at a time
%            default: 65536; if memory is not an issue, use numel(I)
% NumPoints - number of Chebyshev-Lobatto grid points to query within
%                   an interval
%             default: 5
% GlobalMaxScale - true to search for the global scale maximum
%                  false to accept the first maximum found in the interval
%                  default: true
%                      
%
% OUTPUT
% scale - Best scale within the range of scales given
% orientation - Best orientation given the angularOrder
% scaleV - Response at the best scale
% orientationV - Response at the best orientation

% Mark Kittisopikul
% UT Southwestern
% July 14, 2016

% TODO: Review minandmax2 of chebfun2 (separableApprox.minandmax2)

    ip = inputParser;
    ip.addRequired('image',@(x) validateattributes(x,{'numeric'},{'2d'},'image'));
    ip.addParameter('AngularOrder',5, ...
        @(x) validateattributes(x,{'numeric'},{'scalar','integer'}, ...
        'AngularOrder'));
    ip.addParameter('Scales',1:4, ...
        @(x) validateattributes(x,{'numeric'},{'positive','vector'}, ...
        'Scales'));
    ip.addParameter('Mask',true(size(I)), ...
        @(x) validateattributes(x,{'numeric','logical'},{'binary','2d'}, ...
        'Mask'));
    ip.addParameter('TileSize',65536, ...
        @(x) validateattributes(x,{'numeric'},{'scalar','integer'}, ...
        'TileSize'));
    ip.addParameter('NumPoints',5, ...
        @(x) validateattributes(x,{'numeric'},{'scalar','integer'}, ...
        'TileSize'));
    ip.addParameter('GlobalMaxScale',true, ...
        @(x) validateattributes(x,{'numeric','logical'},{'scalar','binary'}, ...
        'GlobalMaxScale'));
    ip.addParameter('NumWorkers',~isempty(gcp('nocreate'))*realmax, ...
        @(x) validateattributes(x,{'numeric'},{'scalar','integer'}, ...
        'NumWorkers'));
    ip.parse(I, varargin{:});
    angularOrder = ip.Results.AngularOrder;
    scales = ip.Results.Scales;
    mask = ip.Results.Mask;
    tileSize = ip.Results.TileSize;
    npts = ip.Results.NumPoints;
    globalMaxScale = ip.Results.GlobalMaxScale;
    ip.delete;

    nAngles = 2*angularOrder+1;
    
    % Initialize outputs
    scale = zeros(size(I));
    orientation = zeros(size(I));
    scaleV = zeros(size(I));
    orientationV = zeros(size(I));
    orientationAtMaxScale = zeros([nAngles size(I)]);
    
    for s = 1:length(scales)-1
        % Scan scales along the Chebyshev-Lobatto grid between pairs of
        % scale values
        F = OrientationSpaceRidgeFilter(1./(2*pi)./chebpts(npts,scales([0 1]+s)),[],angularOrder);
        % Calculate once
        % Or for each scale range if not looking for globalMaxScale
        if(s == 1 || ~globalMaxScale)
            % Number of mask pixels
            nMaskPx = sum(mask(:));
            % Tile to save memory
            nTiles = ceil(nMaskPx/tileSize);
            % Number the mask pixels
            cmask = reshape(cumsum(mask(:)),size(mask));
        end
        % Initialize orientationScaleSpace
        R = F*I;
        % Rearrange to orienation by scale by Y by X
        orientationScaleSpace = permute(R.getArraySpace,[3 4 1 2]);
        orientationScaleSpace = orientationScaleSpace(:,:,mask);
        clear F R A
        for t=1:nTiles
            tileMask = mask & cmask > (t-1)*tileSize & cmask <= t*tileSize;
            
            % Find the best orientation at each discrete scale
            [orientationMax,~,orientationMaxValue] = interpft_extrema(orientationScaleSpace(:,:,(1+(t-1)*tileSize):min(end,t*tileSize)),1,true);
            % Leave option to save the local maxima?
            orientationMax = squeeze(orientationMax(1,:,:));
            orientationMaxValue = squeeze(orientationMaxValue(1,:,:));
            
            % Across the best orientations, find the best scale
            % TODO: Might need to lower the tolerance here
            [scaleMax,~,scaleMaxValue] = interpft_extrema(orientationMaxValue([end:-1:1 2:end-1],:),1,true,1e-9);
            scaleMax = scaleMax(1,:);
            scaleMaxValue = scaleMaxValue(1,:);
            
            % Leave option to save the local maxima?
            [orientationMaxValue,idx] = max(orientationMaxValue);
            orientationMax = orientationMax(sub2ind(size(orientationMax),idx,1:size(orientationMax,2)));
            
            % Map scale from Chebyshev-Lobatto grid on [-1 1] to interval
            % between scale pairs
            interval = scales(s+1) - scales(s);
%             scaleMax = cos(scaleMax)/2*interval+(1/2*interval)+scales(s);
                       
            if(globalMaxScale)
                isGreater = scaleMaxValue(:) > scaleV(tileMask);
                tileMask(tileMask) = isGreater;
                % Save the best orientation
                orientation(tileMask) = orientationMax(isGreater);
                orientationV(tileMask) = orientationMaxValue(isGreater);
                % Interpolate scale to get the orientation
                orientationAtMaxScale(:,tileMask) = orientationSpace.interpolateScale(orientationScaleSpace(:,:,isGreater),scaleMax(isGreater));
                % Save the best scale
                scaleMax = cos(scaleMax)/2*interval+(1/2*interval)+scales(s);
                scale(tileMask) = scaleMax(isGreater);
                scaleV(tileMask) = scaleMaxValue(isGreater);
            else
                % Save the best orientation
                orientation(tileMask) = orientationMax;
                orientationV(tileMask) = orientationMaxValue;
                % Interpolate scale to get the orientation
                orientationAtMaxScale(:,tileMask) = orientationSpace.interpolateScale(orientationScaleSpace,scaleMax);
                % Save the best scale
                scaleMax = cos(scaleMax)/2*interval+(1/2*interval)+scales(s);
                scale(tileMask) = scaleMax;
                scaleV(tileMask) = scaleMaxValue;
                % Mask the values that are at the achieve local max in the interval
                mask(tileMask) = scaleMax == scales(s+1);
            end
        end
    end
    orientationAtMaxScale = permute(orientationAtMaxScale,[2 3 1]);
end