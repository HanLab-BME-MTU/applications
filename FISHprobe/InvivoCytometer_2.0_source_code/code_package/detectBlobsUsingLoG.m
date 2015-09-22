function [imBlobLocations, varargout ] = detectBlobsUsingLoG( im, meanBlobDiameter, varargin )
% Detects blobs as local maxima of the LoG filtered image
%
% The input intensity image is first filtered with a Laplacian of Gaussian (LoG)
% Filter and then the blobs are detected as local maxima in this
% filtered image.
%
% References:
%   
% Byun, J., M. R. Verardo, et al. (2006). 
% "Automated tool for the detection of cell nuclei in digital microscopic 
% images: application to retinal images." Mol Vis 12: 949-960.
%
% Author: Deepak Roy Chittajallu
%

    p = inputParser;    
    p.addRequired( 'im', @(x) (isnumeric(x) && ismember( ndims(x), [2,3] )) );
    p.addRequired( 'meanBlobDiameter', @(x) isscalar(x) );
    p.parse( im, meanBlobDiameter );
    p.addParamValue( 'minBlobDistance', [], @(x) isscalar(x) );
    p.addParamValue( 'spacing', ones( 1, ndims(im) ), @(x) (isnumeric(x) && numel(x) == ndims(im)) );
    p.addParamValue( 'flagBrightBlobs', true, @(x) (isscalar(x) && islogical(x)) );    
    p.addParamValue( 'debugMode', false, @(x) (isscalar(x) && islogical(x)) );
    p.parse( im, meanBlobDiameter, varargin{:} ); 
    
    minBlobDistance = p.Results.minBlobDistance;
    spacing = p.Results.spacing;
    flagBrightBlobs = p.Results.flagBrightBlobs;
    flagDebugMode = p.Results.debugMode;
    
    % Apply Laplacian of Gaussian (LoG) filter
    sigma = 0.5 * meanBlobDiameter / sqrt( ndims(im) );
    imLoG = filterLoGND( im, sigma, 'spacing', spacing, 'UseNormalizedDerivatives', true );
    if flagBrightBlobs
        imLoG = -1 * imLoG;
    end

    % locate local intensity maxima in gaussian blurred image
    if isempty(minBlobDistance)
        MaximaSuppressionSize = 3;
    else
        MaximaSuppressionSize = round( minBlobDistance ./ spacing );    
        evenind = (mod( MaximaSuppressionSize, 2 ) == 0);
        MaximaSuppressionSize( evenind ) = MaximaSuppressionSize( evenind ) + 1;    
        MaximaSuppressionSize( MaximaSuppressionSize < 3 ) = 3; 
    end
    
    switch ndims( im ) 
        
        case 2 
            
            imLocalMax = locmax2d(imLoG, MaximaSuppressionSize, 1);            
    
            if flagDebugMode

                [cy, cx] = find( imLocalMax > 0 );                
                theta = 0:0.1:(2*pi+0.1);
                cx = cx(:,ones(size(theta)));
                cy = cy(:,ones(size(theta)));
                cellRadius = (0.5 * meanBlobDiameter);
                rad = cellRadius * ones( size(cx) );
                theta = theta(ones(size(cx,1),1),:);                
                X = cx + cos(theta).* rad;
                Y = cy + sin(theta).* rad;
                
                imseriesmaskshow( imLoG, imdilate( imLocalMax, strel('disk',3) ) ); 
                set( gcf, 'Name', 'Local Maxima Overlayed on the response of LoG Filter' );
                hold on;               
                    line(X', Y', 'Color', 'g', 'LineWidth', 1.5);                
                hold off;

                imseriesmaskshow( im, imdilate( imLocalMax, strel('disk',3) ) ); 
                set( gcf, 'Name', 'Seed Points Overlayed on Input Image' );
                hold on;               
                    line(X', Y', 'Color', 'g', 'LineWidth', 1.5);                
                hold off;
                
            end
            
        case 3
            
            imLocalMax = locmax3d(imLoG, MaximaSuppressionSize, 'ClearBorder', false);           

            if flagDebugMode
                imseriesmaskshow( imLoG, imdilate(imLocalMax, ones(3,3,3)) ); 
                set( gcf, 'Name', 'Local Maxima overlayed on the response of LoG Filter' );
                
                imseriesmaskshow( im, imdilate(imLocalMax, ones(3,3,3)) ); 
                set( gcf, 'Name', 'Seed Points Overlayed on Input Image' );
            end
            
    end

    % detect local intensity maxima as cell seed points
    imBlobLocations = imLocalMax;

    % detect response if requested
    if nargout > 1 
        varargout{1} = imLoG;
    end
    
end