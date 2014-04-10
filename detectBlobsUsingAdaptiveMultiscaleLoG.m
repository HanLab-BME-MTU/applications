function [ imBlobSeedPoints, varargout ] = detectBlobsUsingAdaptiveMultiscaleLoG( im, foregroundDistanceMap, varargin )
% Detects cell seed points as local maxima of the multiscale LoG filtered
% image
%
% [ imBlobSeedPoints ] = detectBlobsUsingAdaptiveMultiscaleLoG( im, foregroundDistanceMap, varargin )
% [ imBlobSeedPoints, imMultiscaleLoGResponse ] = detectBlobsUsingAdaptiveMultiscaleLoG( im, foregroundDistanceMap, varargin )
% [ imBlobSeedPoints, imMultiscaleLoGResponse, imBlobSize ] = detectBlobsUsingAdaptiveMultiscaleLoG( im, foregroundDistanceMap, varargin )
% 
%  The maximum scale at each pixel is set using a distance map of the
%  binary mask of cell foreground obtained using a thresholding algorithm.
% 
%  This is claimed to be better than detecting blobs using multiscale LoG
%  because setting maximum sigma/scale for each pixel is based on its distance
%  to the nearest background pixel avoids overblurring.
% 
% References:
% 
% Al-Kofahi, Y., W. Lassoued, et al. (2010). 
% "Improved Automatic Detection and Segmentation of Cell Nuclei in 
% Histopathology Images." IEEE Transactions on Biomedical Engineering, 
% 57(4): 841-852.
% 
% Author: Deepak Roy Chittajallu
% 
% 

    p = inputParser;
    p.CaseSensitive = false;
    p.addRequired( 'im', @(x) (isnumeric(x) && ismember( ndims(x), [2,3] )) );
    p.parse( im );
    
    p.addRequired( 'foregroundDistanceMap', @(x) (isnumeric(x) && ndims(x) == ndims(im) && ~any(size(x) ~= size(im))) );    
    p.addParamValue( 'spacing', ones( 1, ndims(im) ), @(x) (isnumeric(x) && numel(x) == ndims(im)) );
    p.addParamValue( 'blobDiameterRange', [], @(x) (numel(x) == 2) );    
    p.addParamValue( 'flagBrightBlobs', true, @(x) (isscalar(x) && islogical(x)) );
    
    p.addParamValue( 'numLoGScales', 10, @(x) isscalar(x) );
    p.addParamValue( 'minBlobDistance', [], @(x) isscalar(x) );    
    p.addParamValue( 'logResponseCutoff', 0, @(x) isscalar(x) );
    p.addParamValue( 'seedToBackgroundDistanceCutoff', 0, @(x) isscalar(x) );
    
    p.addParamValue( 'debugMode', false, @(x) (isscalar(x) && islogical(x)) );
    p.addParamValue( 'showPlots', false, @(x) (isscalar(x) && islogical(x)) );
    p.addParamValue( 'flagParallelize', false, @(x) (isscalar(x) && islogical(x)) );
    p.parse( im, foregroundDistanceMap, varargin{:} );
    
    imDistMap = foregroundDistanceMap;
    blobDiameterRange = p.Results.blobDiameterRange;
    spacing = p.Results.spacing;
    flagBrightBlobs = p.Results.flagBrightBlobs;
    
    minBlobDistance = p.Results.minBlobDistance;
    numLoGScales = p.Results.numLoGScales;
    logResponseCutoff = p.Results.logResponseCutoff; % weak response cutoff
    seedToBackgroundDistanceCutoff = p.Results.seedToBackgroundDistanceCutoff;
    
    flagDebugMode = p.Results.debugMode;
    flagShowPlots = p.Results.showPlots;
    flagParallelize = p.Results.flagParallelize; % more memory will be consumed    

    % if blobDiameterRange is not specified, use distance map to estimate min and max cell diameter
    if isempty( blobDiameterRange )
        
        if flagDebugMode
            fprintf( '\n\n\tEstimating Blob Diameter Range From Distance Map ...\n\n' );
        end
        
        imDMapMaxima = imregionalmax( imDistMap );    
        blobDiameterRange = 2 * [ min(imDistMap(imDMapMaxima)), ...
                                  max(imDistMap(imDMapMaxima))];
        if flagDebugMode
            fprintf( '\n\tBlob Diameter Range: [%.2f, %.2f]\n', blobDiameterRange(1), blobDiameterRange(2) );
            imseriesmaskshow( imDistMap, imDMapMaxima );
            set( gcf, 'Name', 'Distance Map overlayed with points used to estimate cell diameter range' );  
        end       
        
    end
    
    % Run the LoG filter accross the scale space and record the optimal
    % response and the scale at which the optimal response was found for
    % each pixel. The max sigma of each pixel is adapted using the value
    % of the distance map at that location
    maxDiameterFromDistMap = 2.0 * max( imDistMap(:) );
    blobDiameterRangeAdjusted = sort( blobDiameterRange );

    if blobDiameterRangeAdjusted(2) < maxDiameterFromDistMap
        fprintf('\nWARNING: Data seems to contain objects smaller than the specified maximum diameter.\nSpecified Max diameter - %f, Estimated Max diameter - %f\n', max(blobDiameterRange), maxDiameterFromDistMap);
    end
    
    if blobDiameterRangeAdjusted(2) > maxDiameterFromDistMap
        blobDiameterRangeAdjusted(2) = maxDiameterFromDistMap;
        fprintf('\nNote: Data seems to contain objects smaller than the specified maximum diameter.\nSpecified Max diameter - %f, Estimated Max diameter - %f\n', max(blobDiameterRange), maxDiameterFromDistMap);
    end
    
    sigmaLogRange = log2( (0.5 * sort(blobDiameterRangeAdjusted) / sqrt(ndims(im))) );
    sigmaLogValues = linspace( sigmaLogRange(1), sigmaLogRange(2), numLoGScales);
    
%     sigmaLogStep = sigmaLogValues(2) - sigmaLogValues(1);
%     sigmaLogValues = [sigmaLogRange(1) - sigmaLogStep, sigmaLogValues];
    
    sigmaValues = 2.^sigmaLogValues;
    
    if flagDebugMode 
       fprintf( '\nRunning the scale-adaptive LoG filter at multiple scales on an image of size [ %s ] ...\n', ... 
                sprintf( ' %d ', size(im) ) );  
    end
    
    imSigmaMap = imDistMap / sqrt( ndims(im) );
    
    if flagParallelize
       
        cellLoGResponse = cell(1, numel(sigmaValues));
        
        parfor i = 1:numel( sigmaValues )

            if flagDebugMode
                fprintf( '\n\t%d/%d: Trying sigma value of %.2f for cells of diameter %.2f ... ', ...
                         i, numel( sigmaValues ), sigmaValues(i), ...
                         sigmaValues(i) * 2 * sqrt(ndims(im)) );   
                tic
            end

            cellLoGResponse{i} = filterLoGND( im, sigmaValues(i), ... 
                                              'spacing', spacing, ...
                                              'UseNormalizedDerivatives', true );

            if flagDebugMode
                timeElapsed = toc;
                fprintf( 'It took %.2f seconds\n', timeElapsed );           
            end
            
        end

        imMultiscaleLoGResponse = zeros( size(im) );
        pixelScaleMap = ones( size(im) );
        
        for i = 1:numel( sigmaValues )
            
            if i == 1

               imMultiscaleLoGResponse = cellLoGResponse{i};
               pixelScaleMap = ones( size( im ) );

            else            

                imScaleMask = imSigmaMap >= sigmaValues(i);

                if any( imScaleMask(:) )

                    imBetterMask =  imScaleMask & cellLoGResponse{i} < imMultiscaleLoGResponse;
                    imMultiscaleLoGResponse( imBetterMask ) = cellLoGResponse{i}( imBetterMask );
                    pixelScaleMap( imBetterMask ) = i;

                else

                    break;

                end

            end
            
        end
        
    else
        
        for i = 1:numel( sigmaValues )

            if flagDebugMode
                fprintf( '\n\t%d/%d: Trying sigma value of %.2f for cells of diameter %.2f ... ', ...
                         i, numel( sigmaValues ), sigmaValues(i), ...
                         sigmaValues(i) * 2 * sqrt(ndims(im)) );   
                tic
            end

            [ imCurLoGResponse ] = filterLoGND( im, sigmaValues(i), ... 
                                                'spacing', spacing, ...
                                                'UseNormalizedDerivatives', true );

            if flagDebugMode
                timeElapsed = toc;
                fprintf( 'It took %.2f seconds\n', timeElapsed );           
            end

            if i == 1

               imMultiscaleLoGResponse = imCurLoGResponse;
               pixelScaleMap = ones( size( im ) );

            else            

                imScaleMask = imSigmaMap >= sigmaValues(i);

                if any( imScaleMask(:) )

                    imBetterMask =  imScaleMask & imCurLoGResponse < imMultiscaleLoGResponse;
                    imMultiscaleLoGResponse( imBetterMask ) = imCurLoGResponse( imBetterMask );
                    pixelScaleMap( imBetterMask ) = i;

                else

                    break;

                end

            end        

        end
        
    end
    
    if flagBrightBlobs
        imMultiscaleLoGResponse = -1 * imMultiscaleLoGResponse;
    end
    
    % locate local intensity maxima in gaussian blurred image   
    if isempty(minBlobDistance)
        MaximaSuppressionSize = round(  min( blobDiameterRange ) ./ spacing );
    else
        MaximaSuppressionSize = round( minBlobDistance ./ spacing );    
    end

    evenind = (mod( MaximaSuppressionSize, 2 ) == 0);
    MaximaSuppressionSize( evenind ) = MaximaSuppressionSize( evenind ) + 1;    
    MaximaSuppressionSize( MaximaSuppressionSize < 3 ) = 3; 
    
    switch ndims( im ) 
        
        case 2 
            
            imLocalMax = locmax2d(imMultiscaleLoGResponse, MaximaSuppressionSize, 1);            
            
        case 3
            
            imLocalMax = locmax3d( imMultiscaleLoGResponse, MaximaSuppressionSize, ...
                                   'ClearBorder', false);          
            
    end
    
    imLocalMax( imLocalMax > 0 & (imMultiscaleLoGResponse < logResponseCutoff) ) = 0;
    imLocalMax( imLocalMax > 0 & (imDistMap < seedToBackgroundDistanceCutoff) ) = 0;
    
    if flagShowPlots
        
        switch ndims( im ) 

            case 2
                
                seedInd = find( imLocalMax > 0 );
                [cy, cx] = ind2sub( size(imLocalMax), seedInd );                
                theta = 0:0.1:(2*pi+0.1);

                cx = cx(:,ones(size(theta)));
                cy = cy(:,ones(size(theta)));
                rad = 2.0 * (sigmaValues(pixelScaleMap(seedInd)))' * sqrt(ndims(im));
                rad = rad(:,ones(size(theta)));

                theta = theta(ones(size(cx,1),1),:);                
                X = cx + cos(theta).* rad;
                Y = cy + sin(theta).* rad;

                imseriesmaskshow( imMultiscaleLoGResponse, imdilate( imLocalMax, strel('disk',3) ) ); 
                set( gcf, 'Name', 'Local Maxima Overlayed on the response of Adaptive Multiscale LoG Filter' );
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

            case 3

                seedPixelInd = find( imLocalMax > 0 );

                seedPos = cell(1, 3);
                [seedPos{:}] = ind2sub( size(imLocalMax), seedPixelInd );
                seedPos = cat( 2, seedPos{[2,1]}, seedPos{3:end} );

                kd = KDTreeSearcher( seedPos * diag(spacing) );

                pixelPos = cell(1, 3);
                [pixelPos{:}] = ind2sub( size(imLocalMax), (1:numel(imLocalMax))' );
                pixelPos = cat( 2, pixelPos{[2,1]}, pixelPos{3:end} );     

                [closestSeedInd,distanceToSeed] = kd.knnsearch(pixelPos * diag(spacing));
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