function [ imBlobSeedPoints, varargout ] = detectBlobsUsingMultiscaleLoG( im, blobDiameterRange, varargin )
% Detects blobs as local maxima of the Multiscale LoG filtered image
% 
%   [ imBlobLocations ] = detectBlobsUsingMultiscaleLoG( im, blobDiameterRange, varargin )
%   [ imBlobLocations, imMultiscaleLoGResponse ] = detectBlobsUsingMultiscaleLoG( im, blobDiameterRange, varargin )
%   [ imBlobLocations, imMultiscaleLoGResponse, imBlobSize ] = detectBlobsUsingMultiscaleLoG( im, blobDiameterRange, varargin )
% 
% The input intensity image is first filtered with a Laplacian of Gaussian
% (LoG) Filter accross multiple scales/sigmas.
% The seed points are then detected as local maxima in the scale space.
%
% References:
%   
% Byun, J., M. R. Verardo, et al. (2006). 
% "Automated tool for the detection of cell nuclei in digital microscopic 
% images: application to retinal images." Mol Vis 12: 949-960.
%
% Author: Deepak Roy Chittajallu
%
%

    p = inputParser;    
    p.addRequired( 'im', @(x) (isnumeric(x) && ismember( ndims(x), [2,3] )) );
    p.addRequired( 'blobDiameterRange', @(x) (numel(x) == 2) );    
    p.addParamValue( 'spacing', ones( 1, ndims(im) ), @(x) (isnumeric(x) && numel(x) == ndims(im)) );
    p.addParamValue( 'flagBrightBlobs', true, @(x) (isscalar(x) && islogical(x)) );
    p.addParamValue( 'minBlobDistance', [], @(x) isscalar(x) );
    p.addParamValue( 'numLoGScales', 15, @(x) isscalar(x) );
    p.addParamValue( 'logResponseCutoff', eps, @isscalar );
    p.addParamValue( 'roiMask', [], @(x) ( (isnumeric(x) || islogical(x)) && ndims(x) == ndims(im) && all(size(x) == size(im)) ) )    
    p.addParamValue( 'debugMode', false, @(x) (isscalar(x) && islogical(x)) );
    p.parse( im, blobDiameterRange, varargin{:} ); 

    PARAMETERS = p.Results;
    
    % Compute LoG at a series of sigmas
    sigmaLogRange = log2( (0.5 * sort(blobDiameterRange) / sqrt(ndims(im))) );
    sigmaValues = 2.^linspace( sigmaLogRange(1), sigmaLogRange(2), PARAMETERS.numLoGScales);

    [ imMultiscaleLoGResponse, pixelScaleMap ] = filterMultiscaleLoGND( im, sigmaValues, ...
                                                                        'spacing', PARAMETERS.spacing, ...
                                                                        'debugMode', PARAMETERS.debugMode);

    if PARAMETERS.flagBrightBlobs
        imMultiscaleLoGResponse = -1 * imMultiscaleLoGResponse;
    end
    
    % locate local intensity maxima in gaussian blurred image
    if isempty(PARAMETERS.minBlobDistance)
        MaximaSuppressionSize = round(  min( blobDiameterRange ) ./ PARAMETERS.spacing );
    else
        MaximaSuppressionSize = round( PARAMETERS.minBlobDistance ./ PARAMETERS.spacing );    
    end
    
    evenind = (mod( MaximaSuppressionSize, 2 ) == 0);
    MaximaSuppressionSize( evenind ) = MaximaSuppressionSize( evenind ) + 1;    
    MaximaSuppressionSize( MaximaSuppressionSize < 3 ) = 3; 
    
    switch ndims( im ) 
        
        case 2 
            
            imLocalMax = locmax2d(imMultiscaleLoGResponse, MaximaSuppressionSize, 1);            
            
            % prune extrema with weak contrast
            imLocalMax(imLocalMax < PARAMETERS.logResponseCutoff) = 0;

            % prune extrema with weak contrast
            if ~isempty(PARAMETERS.roiMask)
                imLocalMax(~PARAMETERS.roiMask) = 0;
            end
            
            if PARAMETERS.debugMode
                
                seedInd = find( imLocalMax > 0 );
                [cy, cx] = ind2sub( size(imLocalMax), seedInd );                
                theta = 0:0.1:(2*pi+0.1);
                
                cx = cx(:,ones(size(theta)));
                cy = cy(:,ones(size(theta)));
                rad = (sigmaValues(pixelScaleMap(seedInd)))' * sqrt(ndims(im));
                rad = rad(:,ones(size(theta)));
                
                theta = theta(ones(size(cx,1),1),:);                
                X = cx + cos(theta).* rad;
                Y = cy + sin(theta).* rad;

                imseriesmaskshow( imMultiscaleLoGResponse, imdilate( imLocalMax, strel('disk',3) ) ); 
                set( gcf, 'Name', 'Local Maxima Overlayed on the response of Multiscale LoG Filter' );
                hold on;               
                    line(X', Y', 'Color', 'g', 'LineWidth', 1.5);                
                hold off;
                
                imseriesmaskshow( im, imdilate( imLocalMax, strel('disk',3) ) ); 
                set( gcf, 'Name', 'Seed Points Overlayed on Input Image' );
                hold on;               
                    line(X', Y', 'Color', 'g', 'LineWidth', 1.5);                
                hold off;
                
                figure, histogram( pixelScaleMap(seedInd) );
                title( 'Histogram of the cell scales (diameter) found in the image' );
                
            end
            
        case 3
            
            imLocalMax = locmax3d(imMultiscaleLoGResponse, MaximaSuppressionSize, ...
                                  'ClearBorder', false);           
            
            % prune extrema with weak contrast
            imLocalMax(imLocalMax < PARAMETERS.logResponseCutoff) = 0;

            % prune extrema with weak contrast
            if ~isempty(PARAMETERS.roiMask)
                imLocalMax(~PARAMETERS.roiMask) = 0;
            end
                              
            if PARAMETERS.debugMode
                
                seedPixelInd = find( imLocalMax > 0 );
                
                seedPos = cell(1, 3);
                [seedPos{:}] = ind2sub( size(imLocalMax), seedPixelInd );
                seedPos = cat( 2, seedPos{[2,1]}, seedPos{3:end} );
                
                kd = KDTreeSearcher( seedPos * diag(PARAMETERS.spacing) );
                
                pixelPos = cell(1, 3);
                [pixelPos{:}] = ind2sub( size(imLocalMax), (1:numel(imLocalMax))' );
                pixelPos = cat( 2, pixelPos{[2,1]}, pixelPos{3:end} );     

                [closestSeedInd, distanceToSeed] = kd.knnsearch(pixelPos * diag(PARAMETERS.spacing));
                closestSeedPixelInd = seedPixelInd( closestSeedInd );
                
                imBlobMask = zeros( size(imLocalMax) );
                imBlobRadius = sigmaValues( pixelScaleMap ) * sqrt(ndims(im));
                flagIsPixelInSeedVicinity = distanceToSeed <= imBlobRadius(closestSeedPixelInd);
                imBlobMask( flagIsPixelInSeedVicinity ) = closestSeedInd( flagIsPixelInSeedVicinity );
                
                seedMask = imdilate(imLocalMax, ones(3,3,3));
                imseriesmaskshow( imMultiscaleLoGResponse, {seedMask, imBlobMask}); 
                set( gcf, 'Name', 'Local Maxima Overlayed on the response of Adaptive Multiscale LoG Filter' );
                
                imseriesmaskshow( im, {seedMask, imBlobMask}); 
                set( gcf, 'Name', 'Seed Points Overlayed on Input Image' );
                
            end
            
    end

    % detect local intensity maxima as cell seed points
    imBlobSeedPoints = imLocalMax;   
    if nargout > 1
        varargout{1} = imMultiscaleLoGResponse;
    end 
    
    if nargout > 2
        % return estimated scale or radius for each detected blob
        varargout{2} = sigmaValues( pixelScaleMap ) * sqrt(ndims(im)); 
    end
    
end