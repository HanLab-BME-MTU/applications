function [ imCellSeedPoints, varargout ] = detect_cell_seeds_adaptive_multiscale_LoG( im, imCellForegroundMask, varargin )
% Detects cell seed points as local maxima of the multiscale LoG filtered
% image
%
% [ imCellSeedPoints ] = detect_cell_seeds_adaptive_multiscale_LoG( im, imCellForegroundMask, varargin )
% [ imCellSeedPoints, imMultiscaleLoGResponse ] = detect_cell_seeds_adaptive_multiscale_LoG( im, imCellForegroundMask, varargin )
%
%  The maximum scale at each pixel is set using a distance map of the
%  binary mask of cell foreground obtained using a thresholding algorithm
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
    
    p.addRequired( 'imCellForegroundMask', @(x) (isnumeric(x) && ndims(x) == ndims(im) && ~any(size(x) ~= size(im))) );
    p.addParamValue( 'cellDiameterRange', [], @(x) (numel(x) == 2) );
    p.addParamValue( 'spacing', ones( 1, ndims(im) ), @(x) (isnumeric(x) && numel(x) == ndims(im)) );
    p.addParamValue( 'numLoGScales', 10, @(x) isscalar(x) );
    p.addParamValue( 'debugMode', false, @(x) (isscalar(x) && islogical(x)) );
    p.parse( im, imCellForegroundMask, varargin{:} );
    
    cellDiameterRange = p.Results.cellDiameterRange;
    spacing = p.Results.spacing;
    numLoGScales = p.Results.numLoGScales;
    flagDebugMode = p.Results.debugMode;
    
    % compute distance map
    if flagDebugMode
        fprintf( '\n\n\tComputing Distance Map of Binary Foreground Mask ...\n\n' );
    end
    
    imDistMap = bwdistsc( ~imCellForegroundMask, spacing );
    
    if flagDebugMode
        imseriesmaskshow( imDistMap, imCellForegroundMask );
        set( gcf, 'Name', 'Distance Map' );    
    end

    % if not specified, use distance map to estimate min and max cell diameter
    if isempty( cellDiameterRange )
        
        if flagDebugMode
            fprintf( '\n\n\tEstimating Cell Diameter Range From Distance Map ...\n\n' );
        end
        
        imDMapMaxima = imregionalmax( imDistMap );    
        cellDiameterRange = 2 * [ min(imDistMap(imDMapMaxima)), ...
                                  max(imDistMap(imDMapMaxima))];
        if flagDebugMode
            fprintf( '\n\tCell Diameter Range: [%.2f, %.2f]\n', cellDiameterRange(1), cellDiameterRange(2) );
            imseriesmaskshow( imDistMap, imDMapMaxima );
            set( gcf, 'Name', 'Distance Map overlayed with points used to estimate cell diameter range' );  
        end       
        
    end
    
    % Run the LoG filter accross the scale space and record the optimal
    % response and the scale at which the optimal response was found for
    % each pixel. The max sigma of each pixel is adapted using the value
    % of the distance map at that location
    
    sigmaLogRange = log2( (0.5 * sort(cellDiameterRange) / sqrt(ndims(im))) );
    sigmaValues = 2.^linspace( sigmaLogRange(1), sigmaLogRange(2), numLoGScales);

    if flagDebugMode 
       fprintf( '\nRunning LoG filter at multiple scales on an image of size [ %s ] ...\n', ... 
                sprintf( ' %d ', size(im) ) );  
    end
    
    imSigmaMap = imDistMap / sqrt( ndims(im) );
    
    for i = 1:numel( sigmaValues )
        
        if flagDebugMode
            fprintf( '\n\t%d/%d: Trying sigma value of %.2f ... ', i, numel( sigmaValues ), sigmaValues(i) );   
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
    
    imMultiscaleLoGResponse = -1 * imMultiscaleLoGResponse;
    
    % locate local intensity maxima in gaussian blurred image
    minCellDiameterImsp = min( cellDiameterRange ) ./ spacing;
    MaximaSuppressionSize = round( minCellDiameterImsp );
    evenind = (mod( MaximaSuppressionSize, 2 ) == 0);
    MaximaSuppressionSize( evenind ) = MaximaSuppressionSize( evenind ) + 1;    

    switch ndims( im ) 
        
        case 2 
            
            imLocalMax = locmax2d(imMultiscaleLoGResponse, MaximaSuppressionSize, 1);            
            imLocalMax = double(imLocalMax > 0);
            
            if flagDebugMode
                
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
                set( gcf, 'Name', 'Local Maxima Overlayed on the response of Adaptive Multiscale LoG Filter' );
                hold on;               
                    line(X', Y', 'Color', 'g', 'LineWidth', 1.5);                
                hold off;
                
                
                figure, histogram( pixelScaleMap(seedInd) );
                title( 'Histogram of the cell scales (diameter) found in the image' );
                
            end
            
        case 3
            
            imLocalMax = locmax3d(imMultiscaleLoGResponse, MaximaSuppressionSize);           
            imLocalMax = double(imLocalMax > 0);
            
            if flagDebugMode
                imseriesmaskshow( imMultiscaleLoGResponse, imdilate(imLocalMax, ones(3,3,3)) ); 
                set( gcf, 'Name', 'Local Maxima Overlayed on the response of Adaptive Multiscale LoG Filter' );
            end
            
    end

    % detect local intensity maxima as cell seed points
    imCellSeedPoints = imLocalMax;   
    if nargout > 1
        varargout{1} = imMultiscaleLoGResponse;
    end
    
end