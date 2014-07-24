function [ imCellSeedPoints, varargout ] = detect_cell_seeds_radial_symmetry( im, cellRadiusRange, varargin )
% Detects cell seed points using radial symmetry analysis
%
% References:
% 
% Loy, G. and A. Zelinsky (2003). 
% "Fast radial symmetry for detecting points of interest." 
% IEEE Transactions on Pattern Analysis and Machine Intelligence, 
% 25(8): 959-973.
% 
% Example: Radial symmetry analysis for a single circular blob
% 
%     imsize = [60,60];
%     [y,x] = meshgrid(1:imsize(1),1:imsize(2));
%     xc = round(0.5 * imsize(2));
%     yc = round(0.5 * imsize(1));
%     rad = 0.4 * min(imsize);
%     im = 255 * double((x-xc).^2 + (y-yc).^2 - rad^2 <= 0);
%     mask = double((x-xc).^2 + (y-yc).^2 - rad^2 <= 0.1);
% 
%     imCellSeedPoints = detect_cell_seeds_radial_symmetry( im, [rad, rad], 'roiMask', mask, 'numRadiusSamples', 1 );
%     imseriesmaskshow( im, imdilate(imCellSeedPoints, ones(3,3)) );  
% 
% Example: Radial symmetry analysis for a single elliptical blob
% 
%     imsize = [60,60];
%     [x,y] = meshgrid(1:imsize(2),1:imsize(1));
%     xc = round(0.5 * imsize(2));
%     yc = round(0.5 * imsize(1));
%     majorRad = 0.4 * min(imsize);
%     minorRad = 0.3 * min(imsize);
%     im = 255 * double((x-xc).^2/majorRad^2 + (y-yc).^2/minorRad^2 - 1 <= 0);
%     mask = double((x-xc).^2/majorRad^2 + (y-yc).^2/minorRad^2 - 1 <= 0.1);
% 
%     imCellSeedPoints = detect_cell_seeds_radial_symmetry( im, [minorRad-3, minorRad+3], 'roiMask', mask, 'numRadiusSamples', 4 );
%     imseriesmaskshow( im, imdilate(imCellSeedPoints, ones(3,3)) );  
% 
% Author: Deepak Roy Chittajallu   
% 
%

    p = inputParser;
    p.addRequired( 'im', @(x) (isnumeric(x) && ismember( ndims(x), [2,3] )) );
    p.addRequired( 'cellRadiusRange', @(x) (numel(x) == 2) );    
    p.parse( im, cellRadiusRange );    
    
    p.addOptional( 'roiMask', ones(size(im)), @(x) (ndims(x) == ndims(im) && ~any(size(x) ~= size(im))));    
    p.addParamValue( 'spacing', ones( 1, ndims(im) ), @(x) (isnumeric(x) && numel(x) == ndims(im)) );
    p.addParamValue( 'numRadiusSamples', 10, @(x) isscalar(x) );
    p.addParamValue( 'debugMode', true, @(x) (isscalar(x) && islogical(x)) );
    p.addParamValue( 'flagParallelize', true, @(x) (isscalar(x) && islogical(x)) );
    p.parse( im, cellRadiusRange, varargin{:} ); 

    roiMask = logical( p.Results.roiMask > 0 );
    spacing = p.Results.spacing;
    numRadiusSamples = p.Results.numRadiusSamples;
    flagDebugMode = p.Results.debugMode;    
    flagParallelize = p.Results.flagParallelize;
    
    if flagParallelize
        flagPoolOpenedAlready = matlabpool( 'size' ) > 0;        
        if ~flagPoolOpenedAlready 
            matlabpool open;
        end            
    end
    
    % compute image gradient 
    imSmooth = filterGaussND( im, 1.0, 'spacing', spacing );
    imGradVec = cell(1,ndims(im));       
    switch ndims(im)
        
        case 2
            
            sobelKernelY = -1 * fspecial( 'sobel' );
            imGradVec{1} = imfilter( imSmooth, sobelKernelY ) / (2*spacing(1));
            imGradVec{2} = imfilter( imSmooth, sobelKernelY' ) / (2*spacing(2));
            
        case 3

            sobelKernelZ = zeros(3,3,3);
            sobelKernelZ(:,:,1) = [ -1 -3 -1; -3 -6 -3; -1 -3 -1 ];
            sobelKernelZ(:,:,3) = -1 * sobelKernelZ(:,:,1);
            
            imGradVec{1} = imfilter( imSmooth, permute(sobelKernelZ, [3,1,2]) ) / (2*spacing(1));
            imGradVec{2} = imfilter( imSmooth, permute(sobelKernelZ, [2,3,1]) ) / (2*spacing(2));
            imGradVec{3} = imfilter( imSmooth, sobelKernelZ ) / (2*spacing(3));
            
    end
    
    imGradMag = zeros(size(im));
    for i = 1:ndims(im)
        imGradMag = imGradMag + (imGradVec{i}).^2;
    end
    imGradMag = sqrt( imGradMag );
    
    % perform radial summetry analysis
    radiusValueVec = linspace( cellRadiusRange(1), cellRadiusRange(2), numRadiusSamples);
    
    if flagParallelize

        imRSTransform = zeros(size(im));
        parfor rid = 1:numel(radiusValueVec)            
            imRSTransform = imRSTransform + ComputeRSTransform( imGradMag, imGradVec, roiMask, radiusValueVec(rid), spacing );            
        end
        imRSTransform = imRSTransform / numel(radiusValueVec);
        
    else
        
        if flagDebugMode && ndims(im) == 2
            
            imRSTransform = [];
            for rid = 1:numel(radiusValueVec)          
                imRSTransform = cat( 3, imRSTransform, ComputeRSTransform( imGradMag, imGradVec, roiMask, radiusValueVec(rid), spacing ) );            
            end
            
            imseriesmaskshow( imRSTransform, repmat(roiMask, [1,1,numRadiusSamples]) );
            set( gcf, 'Name', 'Radial Symmetry Tranform at all radial samples' );
            
            imRSTransform = mean(imRSTransform, 3);
            
        else

            imRSTransform = zeros(size(im));
            for rid = 1:numel(radiusValueVec)          
                imRSTransform = imRSTransform + ComputeRSTransform( imGradMag, imGradVec, roiMask, radiusValueVec(rid), spacing );            
            end
            imRSTransform = imRSTransform / numel(radiusValueVec);
            
        end
        

    end
    
    % compute local maxima and do non-maximal suppresson
    minCellDiameterImsp = min( 0.5 * cellRadiusRange ) ./ spacing;
    MaximaSuppressionSize = round( minCellDiameterImsp );
    evenind = (mod( MaximaSuppressionSize, 2 ) == 0);
    MaximaSuppressionSize( evenind ) = MaximaSuppressionSize( evenind ) + 1;    
    
    switch ndims(im)
        
        case 2
            
            imLocalMax = locmax2d(imRSTransform, MaximaSuppressionSize, 1);            
            imLocalMax = double(imLocalMax > 0 & roiMask > 0);
            
        case 3
        
            imLocalMax = locmax3d(imRSTransform, MaximaSuppressionSize);            
            imLocalMax = double(imLocalMax > 0 & roiMask > 0);

    end

    if flagDebugMode
        imseriesmaskshow( imRSTransform, { imdilate(imLocalMax, ones(3*ones(1,ndims(im)))), roiMask } ); 
        set( gcf, 'Name', 'Local Maxima Overlayed on radial symmetry transform' );
    end
    
    % detect local intensity maxima as cell seed points
    imCellSeedPoints = imLocalMax;   
    if nargout > 1
        varargout{1} = imRSTransform;
    end 
   
    if flagParallelize && ~flagPoolOpenedAlready
        matlabpool close;
    end    
    
end

function [imRSTransform] = ComputeRSTransform( imGradMag, imGradVec, roiMask, voteRadius, spacing )

    imdims = ndims(imGradMag);
    imsize = size(imGradMag);
    
    GradRange = [ min(imGradMag(roiMask)), max(imGradMag(roiMask)) ];
    gradThresh = GradRange(1) + 0.1 * (GradRange(2) - GradRange(1));

    roiMask = roiMask & (imGradMag > gradThresh);
    pixind = find( roiMask > 0 );
    pixsubind = cell(1,imdims);    
    [pixsubind{:}] = ind2sub( size(roiMask), pixind );    
    
    imOrientationVote = zeros(imsize);
    imMagnitudeVote = zeros(imsize);
    
    for i = 1:numel(pixind)        
        
        pid = pixind(i);
        gvec = zeros(1,imdims);
        ptCur = zeros(1,imdims);
        for j = 1:imdims
           gvec(j) = imGradVec{j}(pid) ./ imGradMag(pid); 
           ptCur(j) = pixsubind{j}(i);
        end

        ptPlusCenter  = round( ptCur + voteRadius * gvec );
        if ~( any( ptPlusCenter < 1 ) || any(ptPlusCenter > imsize) )
            ptCurCell = num2cell( ptPlusCenter );
            curPixInd = sub2ind( imsize, ptCurCell{:} );
            imOrientationVote( curPixInd ) = imOrientationVote( curPixInd ) + 1;
            imMagnitudeVote( curPixInd ) = imMagnitudeVote( curPixInd ) + imGradMag(curPixInd);
        end

        ptMinusCenter = round( ptCur - voteRadius * gvec );            
        if ~( any( ptMinusCenter < 1 ) || any(ptMinusCenter > imsize) )
            ptCurCell = num2cell( ptMinusCenter );
            curPixInd = sub2ind( imsize, ptCurCell{:} );
            imOrientationVote( curPixInd ) = imOrientationVote( curPixInd ) + 1;
            imMagnitudeVote( curPixInd ) = imMagnitudeVote( curPixInd ) + imGradMag(curPixInd);
        end

    end    

    imOrientationVote = 0.1 * imOrientationVote; % normalization
    imMagnitudeVote = 0.1 * imMagnitudeVote; % normalization                
    imRSTransform = filterGaussND( imOrientationVote.^2, 0.5 * voteRadius, 'spacing', spacing );

end