function [imCellSeedPoints, varargout] = detect_cell_seeds_IntensityMaxima( im, meanCellDiameter, varargin )
% Detects cell seed points as local intensity maxima of gaussian blurred
% image
%
% [imCellSeedPoints] = detect_cell_seeds_IntensityMaxima( im, meanCellDiameter, varargin )
% [imCellSeedPoints, imBlurred] = detect_cell_seeds_IntensityMaxima( im, meanCellDiameter, varargin )
% 
    p = inputParser;
    p.addRequired( 'im', @(x) (isnumeric(x) && ismember( ndims(x), [2,3] )) );
    p.parse( im);
    
    imdims = ndims(im);

    p.addRequired( 'meanCellDiameter', @(x) (isnumeric(x) && (isscalar(x) || numel(x) == imdims)) );
    p.addParamValue( 'spacing', ones( 1, ndims(im) ), @(x) (isnumeric(x) && numel(x) == ndims(im)) ); 
    p.addParamValue( 'debugMode', false, @(x) ( islogical(x) && isscalar(x) ) );
    p.parse( im, meanCellDiameter, varargin{:} );
    
    spacing = p.Results.spacing;
    flagDebugMode = p.Results.debugMode;
    meanCellDiameterImsp = meanCellDiameter ./ spacing;
    
    % smooth image before looking for locate maxima
    sigma = meanCellDiameter / 3.5;
    imBlurred = filterGaussND(im, sigma, 'spacing', spacing);

    % locate local intensity maxima in gaussian blurred image
    MaximaSuppressionSize = round( meanCellDiameterImsp );
    evenind = (mod( MaximaSuppressionSize, 2 ) == 0);
    MaximaSuppressionSize( evenind ) = MaximaSuppressionSize( evenind ) + 1;    
    
    switch imdims 
        
        case 2 
            
            imLocalMax = locmax2d(imBlurred, MaximaSuppressionSize);            
            imLocalMax = double(imLocalMax > 0);
            
            if flagDebugMode
                imseriesmaskshow( imBlurred, imdilate( imLocalMax, strel('disk',3) ) ); 
                set( gcf, 'Name', 'Local Maxima overlayed on Blurred Image' );
            end
            
        case 3
            
            imLocalMax = locmax3d(imBlurred, MaximaSuppressionSize);           
            imLocalMax = double(imLocalMax > 0);
            
            if flagDebugMode
                imseriesmaskshow( imBlurred, imdilate(imLocalMax, ones(3,3,3)) ); 
                set( gcf, 'Name', 'Local Maxima overlayed on Blurred Image' );
            end
            
    end

    % detect local intensity maxima as cell seed points
    imCellSeedPoints = imLocalMax;
    if nargout > 1 
        varargout{1} = imBlurred;
    end

end