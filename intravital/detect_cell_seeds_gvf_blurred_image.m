function [imCellSeedPoints] = detect_cell_seeds_gvf_blurred_image( im, minCellDiameter, varargin )

    p = inputParser;
    p.addRequired( 'im', @(x) (isnumeric(x) && ismember( ndims(x), [2,3] )) );
    p.parse( im );  
    
    imdims = ndims(im);

    p.addRequired( 'minCellDiameter', @(x) (isnumeric(x) && (isscalar(x) || numel(x) == imdims)) );
    p.addParamValue( 'maxGVFIterations', 100, @(x) (isnumeric(x) && isscalar(x)) );
    p.addParamValue( 'weightGVFRegularization', 0.15, @(x) (isnumeric(x) && isscalar(x)) );
    p.addParamValue( 'thresholdDivergence', 0.1, @(x) (isnumeric(x) && isscalar(x)) );
    p.addParamValue( 'debugMode', false, @(x) ( islogical(x) && isscalar(x) ) );
    p.parse( im, minCellDiameter, varargin{:} );
    
    flagDebugMode = p.Results.debugMode;
    maxGVFIterations = p.Results.maxGVFIterations;
    weightGVFRegularization = p.Results.weightGVFRegularization;
    thresholdDivergence = p.Results.thresholdDivergence;
    
    switch imdims
        
        case 2

            % compute feature image
            imFeature = filterGaussND(im, minCellDiameter / 3.5);
            
            if flagDebugMode
                figvec = [];
            end
            
            % display intial vector field
            [fx, fy] = gradient( imFeature );
            if flagDebugMode
                imseriesshow( imFeature );
                hold on, quiver(fx,fy), hold off
                set( gcf, 'Name', 'Initial Vector Field' );
                figvec(end+1) = gcf;
            end
            
            % compute diffused gradient vector field
            [vx,vy] = gvfc(imFeature, weightGVFRegularization, maxGVFIterations);            
            vmag = sqrt( vx .* vx + vy .* vy );
            
            if flagDebugMode
                imseriesshow( imFeature );
                hold on, quiver(vx ./ vmag, vy ./ vmag), hold off
                set( gcf, 'Name', 'Diffused Vector Field' );
                figvec(end+1) = gcf;
            end
            
            % compute divergence of the diffused GVF
            imdiv = divergence(vx ./ vmag, vy ./ vmag);

            if flagDebugMode
                imseriesshow( imdiv );                
                set( gcf, 'Name', 'Divergence of Diffused Vector Field' );
                figvec(end+1) = gcf;
            end
            
            % detect local extrema in divergence image
            MaximaSuppressionSize = ceil( minCellDiameter / 1.5 );
            evenind = (mod( MaximaSuppressionSize, 2 ) == 0);
            MaximaSuppressionSize( evenind ) = MaximaSuppressionSize( evenind ) + 1;    
            imLocalMax = locmax2d(abs(imdiv), MaximaSuppressionSize);                

            % detect cell seed points as strong sinks within local extrema of the divergence image
            imCellSeedPoints = (imLocalMax > 0 & imdiv < -thresholdDivergence); 

            if flagDebugMode
                
                imseriesshow( im );  
                set( gcf, 'Name', 'Cell Seed Point Candidates' );
                figvec(end+1) = gcf;                               
                
                for i = 1:numel(figvec)
                    
                    figure( figvec(i) );
                    hold on;
                    [yind,xind] = ind2sub( size(im), find( imLocalMax > 0 ) );
                    plot( xind, yind, 'r+', 'MarkerSize', 10.0 );

                    [yind,xind] = ind2sub( size(im), find( imCellSeedPoints > 0 ) );
                    plot( xind, yind, 'go', 'MarkerSize', 10.0 );
                    hold off;
                    
                end
                
            end
            
        case 3
            
            % compute feature image
            imFeature = filterGaussND(im, minCellDiameter / 3.5);

            % compute diffused gradient vector field
            [vx,vy] = GVFND(gmag, weightGVFRegularization, maxGVFIterations, 'flagDebugMode', flagDebugMode);            
            vmag = sqrt( vx .* vx + vy .* vy + vz .* vz );
            
            % compute divergence of the diffused GVF
            imdiv = divergence(vx ./ vmag, vy ./ vmag, vz ./ vmag);
            
            % detect local extrema in divergence image            
            MaximaSuppressionSize = ceil( minCellDiameter / 1.5 );
            evenind = (mod( MaximaSuppressionSize, 2 ) == 0);
            MaximaSuppressionSize( evenind ) = MaximaSuppressionSize( evenind ) + 1;    
            imLocalMax = locmax3d(abs(imdiv), MaximaSuppressionSize);                
    
            % detect strong sources/sinks as local maxima in divergence image
            imCellSeedPoints = (imLocalMax > 0 & imdiv < -thresholdDivergence); 
    end
    
end