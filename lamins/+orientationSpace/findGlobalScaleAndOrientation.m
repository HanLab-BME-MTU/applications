function [ scale, orientation, scaleV, orientationV ] = findGlobalScaleAndOrientation( I, angularOrder, scales, mask, tileSize, npts)
%FINDGLOBALSCALEANDORIENTATION Finds the scales and orientation using OrientationSpace
%filtering. Scans scale using Chebyshev polynomials
%
% INPUT
% I - image as a double matrix
% angularOrder - see OrientationSpaceFilter
% scales - list of scales as a vector or length at least 2. Adjacent pairs
%          of scales will be used as intervals over which to analyze scale
% mask - initial mask identifying the area of interest
% tileSize - (optional) number of pixels to process at a time
%            default: 65536; if memory is not an issue, use numel(I)
% npts - (optional) number of Chebyshev-Lobatto grid points to query within
%                   an interval
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

    if(nargin < 4 || isempty(mask))
        mask = true(size(I));
    elseif(isscalar(mask))
        initScale = scales(mask);
        F = OrientationSpaceFilter(1/2/pi/initScale,[],angularOrder);
        R = F*I;
        mask = lamins.functions.maskFromSteerable(R);
        clear initScale F R
    end
    if(nargin < 5 || isempty(tileSize))
        tileSize = 65536;
    end
    if(nargin < 6 || isempty(npts))
        npts = 5;
    end
    nAngles = 2*angularOrder+1;
    
    % Initialize outputs
    scale = zeros(size(I));
    orientation = zeros(size(I));
    scaleV = zeros(size(I));
    orientationV = zeros(size(I));
    
    for s = 1:length(scales)-1
        % Scan scales along the Chebyshev-Lobatto grid between pairs of
        % scale values
        F = OrientationSpaceRidgeFilter(1./(2*pi)./chebpts(npts,scales([0 1]+s)),[],angularOrder);
        % Number of mask pixels
        nMaskPx = sum(mask(:));
        % Initialize orientationScaleSpace
        orientationScaleSpace = zeros(nAngles,npts,nMaskPx);
        for f = 1:npts
            R = F(f)*I;
            A = permute(R.a,[3 4 1 2]);
            % Collapse spatial dimensions using mask
            orientationScaleSpace(:,f,:) = A(:,1,mask);
        end
        clear F R A
        % Tile to save memory
        nTiles = ceil(nMaskPx/tileSize);
        % Number the mask pixels
        cmask = reshape(cumsum(mask(:)),size(mask));
        for t=1:nTiles
%             orientationScaleSubSpace = orientationScaleSpace(:,:,(1+(t-1)*tileSize):min(end,t*tileSize));
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
            orientation(tileMask) = orientationMax;
            orientationV(tileMask) = orientationMaxValue;
            
            % Previously, we rescanned the whole matrix to find the best
            % scale
%             [scaleMax,~,scaleMaxValue] = interpft_extrema(orientationScaleSubSpace(:,[end:-1:1 2:end-1],:),2,true);
%             scaleMax = squeeze(scaleMax(:,1,:));
%             scaleMaxValue = squeeze(scaleMaxValue(:,1,:));
%             [scaleMaxValue,idx] = max(scaleMaxValue);
%             scaleMax = scaleMax(sub2ind(size(scaleMax),idx,1:size(scaleMax,2)));

            % Map scale from Chebyshev-Lobatto grid on [-1 1] to interval
            % between scale pairs
            interval = scales(s+1) - scales(s);
            scaleMax = cos(scaleMax)/2*interval+(1/2*interval)+scales(s);
            
            % Save the best scale
            scale(tileMask) = scaleMax;
            scaleV(tileMask) = scaleMaxValue;
            
            % Mask the values that are at the best scale in the interval
            mask(tileMask) = scaleMax == scales(s+1);
        end
    end
end


