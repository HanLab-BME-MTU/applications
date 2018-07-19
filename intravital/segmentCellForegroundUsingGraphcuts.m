function [ imCellForegroundMask ] = segmentCellForegroundUsingGraphCuts( im, imForegroundSeedRegion, imBackgroundSeedRegion, varargin )

    p = inputParser;
    p.addRequired( 'im', @(x) (isnumeric(x) && ismember( ndims(x), [2,3] )) );
    p.parse( im );
    p.addRequired( 'imForegroundSeedRegion', @(x) (isnumeric(x) && ndims(x) == ndims(im) && ~any(size(x) ~= size(im))) );
    p.addRequired( 'imBackgroundSeedRegion', @(x) (isnumeric(x) && ndims(x) == ndims(im) && ~any(size(x) ~= size(im))) )
    p.addParamValue( 'spacing', ones( 1, ndims(im) ), @(x) (isnumeric(x) && numel(x) == ndims(im)) );
    p.addParamValue( 'regularizationWeight', 10.0, @(x) (isscalar(x)) );
    p.addParamValue( 'debugMode', false, @(x) (isscalar(x) && islogical(x)) );
    p.parse( im, varargin{:} );

    imForegroundSeedRegion = p.Results.imForegreoundSeedRegion;
    imBackgroundSeedRegion = p.Results.imBackgroundSeedRegion;
    
    spacing = p.Results.spacing;
    flagDebugMode = p.Results.debugMode;    
    regularizationWeight = p.Results.regularizationWeight;
    
    imsize = size( im );
    
    % Build foreground and background pdf
    fgmodel.mean = mean( im(imForegroundSeedRegion > 0) );
    fgmodel.covariance = var( im(imForegroundSeedRegion > 0) );
    fgmodel.gmobj = gmdistribution( fgmodel.mean, fgmodel.covariance );
    
    bgmodel.mean = mean( im(imBackgroundSeedRegion > 0) );
    bgmodel.covariance = var( im(imBackgroundSeedRegion > 0) );
    bgmodel.gmobj = gmdistribution( bgmodel.mean, bgmodel.covariance );

    if flagDebugMode
        DisplayModelPdf( fgmodel, bgmodel );
    end
    
    % Setup unary clique potentials
    
        % Background
        FirstOrderCliquePotential_Bgnd = reshape( -log( pdf( bgmodel.gmobj , im(:) ) ) , imsize );

            % Fix Inf vals before normalizing
            maxBgndVal = max( FirstOrderCliquePotential_Bgnd( ~isinf( FirstOrderCliquePotential_Bgnd ) ) );
            FirstOrderCliquePotential_Bgnd( isinf( FirstOrderCliquePotential_Bgnd ) ) = maxBgndVal + 1;

            if flagDebugMode
                imseriesshow( FirstOrderCliquePotential_Bgnd );
                title( 'Unary Clique Potential -- Background', 'FontSize' , 12 , 'FontWeight' , 'Bold' );
                set( gcf , 'Name' , 'Unary Clique Potential -- Background' );         
            end
            
    
        % Foreground
        FirstOrderCliquePotential_Fgnd = reshape( -log( pdf( fgmodel.gmobj , im(:) ) ) , imsize );

            % Fix Inf vals before normalizing
            maxFgndVal = max( FirstOrderCliquePotential_Fgnd( ~isinf( FirstOrderCliquePotential_Fgnd ) ) );
            FirstOrderCliquePotential_Fgnd( isinf( FirstOrderCliquePotential_Fgnd ) ) = maxFgndVal + 1;

            if flagDebugMode
                imseriesshow( FirstOrderCliquePotential_Fgnd );
                title( 'Unary Clique Potential -- Foreground', 'FontSize' , 12 , 'FontWeight' , 'Bold' );
                set( gcf , 'Name' , 'Unary Clique Potential -- Foreground' );         
            end
            
    FirstOrderCliquePotential = [ FirstOrderCliquePotential_Bgnd(:) , FirstOrderCliquePotential_Fgnd(:) ];
    
    % Setup pairwise clique potentials
    directionLabels = { 'Y', 'X', 'Z' };
    ndir_enode1 = [];
    ndir_enode2 = [];
    ndir_ecap = [];
        
        % setup pairwise cliques along each dimension and the associate clique potentials        
        img_pixind = (1:numel(im))';
        img_pixsubind = cell(1,ndims(im));
        [img_pixsubind{:}] = ind2sub( imsize, img_pixind );
        
        for i = 1:ndims(im)
            
            % edge cost
            imGaussGrad = filterGaussGradND(im, 1.0, i, 'spacing', spacing );
            sigmaGrad = std( imGaussGrad(~imThreshForegroundMask) );            
            imGradCost = exp( -imGaussGrad.^2 / (2 * sigmaGrad^2) );
            
            % from node indices
            from_node_ind = img_pixind( (img_pixsubind{i} + 1) <= imsize(i) );
            
            % to node indices
            img_neighpixsubind = cell(1,ndims(im));
            for j = 1:ndims(im)
                img_neighpixsubind{j} = img_pixsubind{j}( from_node_ind );
            end
            img_neighpixsubind{i} = img_neighpixsubind{i} + 1;
            
            to_node_ind = sub2ind( imsize, img_neighpixsubind{:} );
            
            % add to list of pairwise cliques
            ndir_enode1 = [ ndir_enode1 ; from_node_ind ];
            ndir_enode2 = [ ndir_enode2 ; to_node_ind ];
            ndir_ecap   = [ ndir_ecap   ; imGradCost(from_node_ind) ];

            if flagDebugMode
                imseriesshow( imGradCost );
                set( gcf , 'Name' , sprintf( 'Second-Order Clique Potentials -- %s-Dir', directionLabels{i} ) );         
            end
            
        end
        
        SecondOrderCliques = [ ndir_enode1 , ndir_enode2 ];
        SecondOrderCliquePotential = zeros( 2 , 2 , numel( ndir_enode1 ) );
        SecondOrderCliquePotential( 1 , 2 , : ) = ndir_ecap;
        SecondOrderCliquePotential( 2 , 1 , : ) = ndir_ecap;            
        SecondOrderCliquePotential = regularizationWeight * SecondOrderCliquePotential;
        
        clear img_pixind;
        clear img_pixsubind;
    
    % Solve Maxflow using K-B maxflow wrapper
    if flagDebugMode
        fprintf( 1 , '\n\nSolving the Maximum Flow Problem using K-B Maxflow wrapper ... '  );
    end

        tic
        [ energy , seg_labels ] = MinimizeBinaryOrder2MRFEnergy_GraphCuts( FirstOrderCliquePotential, SecondOrderCliques, SecondOrderCliquePotential );       
        timeElapsed = toc;

        if flagDebugMode
            fprintf( 1 , 'It took %.2f seconds\n\n', timeElapsed  );
        end
        
        imCellForegroundMask = reshape( seg_labels , imsize );
        imCellForegroundMask = imCellForegroundMask > 0;
    
       % Display
       if flagDebugMode           
           
           imseriesmaskshow( im , { imCellForegroundMask, imThreshForegroundMask } );

           title( 'Cell Foreground Segmentation using MRF' , ...
                  'FontSize' , 12 , 'FontWeight' , 'Bold' );    
           set( gcf , 'Name' , 'Cell Foreground Segmentation using MRF: SegResult -- ThreshResult' );
              
       end
    
end