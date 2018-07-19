function [ imLabeledMask ] = segmentCellsInteractivelyUsingTopoGraphcut( imInput, varargin )

    p = inputParser;
    p.addRequired( 'imInput', @(x) (isnumeric(x) && ismember( ndims(x), [2,3] )) );
    p.parse( imInput );
    p.addParamValue( 'regularizationWeight', 1.0, @(x) (isscalar(x)) );
    p.addParamValue( 'spacing', ones( 1, ndims(imInput) ), @(x) (isnumeric(x) && numel(x) == ndims(imInput)) );
    p.addParamValue( 'debugMode', true, @(x) (isscalar(x) && islogical(x)) );
    p.parse( imInput, varargin{:} );
    
    imInput = double(p.Results.imInput);
    spacing = p.Results.spacing;
    flagDebugMode = p.Results.debugMode;  
    regularizationWeight = p.Results.regularizationWeight;
    
    SEED_NEIGH = 5;
    imsize = size(imInput);    
    imAdjusted = filterGaussND( imInput, 5.0, 'spacing', spacing );
    ImageIntensityRange = ComputeImageDynamicRange( imAdjusted, 98.0 );    

    % get topololgy mask
    imTopologyMask = zeros(size(imAdjusted));    
    
    switch ndims(imInput)
       
        case 2
            
            [ fgnd_seed_points , bgnd_seed_points ] = get_fgnd_bgnd_seeds_2d_points( imInput, ImageIntensityRange );

            for i = 1:size(fgnd_seed_points,1)
                imTopologyMask( fgnd_seed_points(i,2), fgnd_seed_points(i,1) ) = 1;
            end

            imTopologyMask = imdilate( imTopologyMask, ones(3,3) );            

        case 3        

            [ fgnd_seed_points , bgnd_seed_points ] = get_fgnd_bgnd_seeds_3d_points( imInput, ImageIntensityRange );

            for i = 1:size(fgnd_seed_points,1)
                imTopologyMask( fgnd_seed_points(i,2), fgnd_seed_points(i,1), fgnd_seed_points(i,3) ) = 1;
            end

            imTopologyMask = imdilate( imTopologyMask, ones(3,3,3) );
            
    end
    
    if flagDebugMode
        
        imseriesmaskshow( imInput, imTopologyMask );
        title( 'Topology Mask', 'FontSize' , 12 , 'FontWeight' , 'Bold' );
        set( gcf , 'Name' , 'Topology Mask' );        

    end
    
    % get seed regions
    switch ndims(imInput)
       
        case 2
            
            [ fgnd_seed_points , bgnd_seed_points ] = get_fgnd_bgnd_seeds_2d_strokes( imInput, ImageIntensityRange );

            % Build foreground and background pdf
            fgmodel = BuildModelPdf( fgnd_seed_points , imAdjusted, SEED_NEIGH );
            bgmodel = BuildModelPdf( bgnd_seed_points , imAdjusted , SEED_NEIGH );    

        case 3        

            [ fgnd_seed_points , bgnd_seed_points ] = get_fgnd_bgnd_seeds_3d_strokes( imInput, ImageIntensityRange );

            % Build foreground and background pdf
            fgmodel = BuildModelPdf( fgnd_seed_points , imAdjusted, SEED_NEIGH );
            bgmodel = BuildModelPdf( bgnd_seed_points , imAdjusted , SEED_NEIGH );    
            
    end

    if flagDebugMode
        
        DisplayModelPdf( fgmodel , bgmodel );

        imseriesmaskshow( imInput, {fgmodel.im_seed_mask,  bgmodel.im_seed_mask} );
        title( 'Seed Overlay', 'FontSize' , 12 , 'FontWeight' , 'Bold' );
        set( gcf , 'Name' , 'Seed Overlay' );        

    end
    
    %% Setup unary clique potentials
    
        % Background
        FirstOrderCliquePotential_Bgnd = reshape( -log( pdf( bgmodel.gmobj , imAdjusted(:) ) ) , imsize );

            % Fix Inf vals before normalizing
            maxBgndVal = max( FirstOrderCliquePotential_Bgnd( ~isinf( FirstOrderCliquePotential_Bgnd ) ) );
            FirstOrderCliquePotential_Bgnd( isinf( FirstOrderCliquePotential_Bgnd ) ) = maxBgndVal + 1;

            if flagDebugMode
                imseriesshow( FirstOrderCliquePotential_Bgnd );
                title( 'Unary Clique Potential -- Background', 'FontSize' , 12 , 'FontWeight' , 'Bold' );
                set( gcf , 'Name' , 'Unary Clique Potential -- Background' );         
            end
    
        % Foreground
        FirstOrderCliquePotential_Fgnd = reshape( -log( pdf( fgmodel.gmobj , imAdjusted(:) ) ) , imsize );

            % Fix Inf vals before normalizing
            maxFgndVal = max( FirstOrderCliquePotential_Fgnd( ~isinf( FirstOrderCliquePotential_Fgnd ) ) );
            FirstOrderCliquePotential_Fgnd( isinf( FirstOrderCliquePotential_Fgnd ) ) = maxFgndVal + 1;

            if flagDebugMode
                imseriesshow( FirstOrderCliquePotential_Fgnd );
                title( 'Unary Clique Potential -- Foreground', 'FontSize' , 12 , 'FontWeight' , 'Bold' );
                set( gcf , 'Name' , 'Unary Clique Potential -- Foreground' );         
            end    
            
    %% Setup pairwise clique potentials
    directionLabels = { 'Y', 'X', 'Z' };
    ndir_enode1 = [];
    ndir_enode2 = [];
    ndir_ecap = [];
    
        % setup pairwise cliques along each dimension and the associate clique potentials        
        SecondOrderCliquePotential = zeros( numel(imInput), 2*ndims(imInput) );
        img_pixind = (1:numel(imInput))';
        img_pixsubind = cell(1,ndims(imInput));
        [img_pixsubind{:}] = ind2sub( imsize, img_pixind );
        
        for i = 1:ndims(imInput)
            
            % edge cost
            imGaussGrad = filterGaussGradND(imInput, 1.0, i, 'spacing', spacing );
            sigmaGrad = std( imGaussGrad(fgmodel.im_seed_mask) );            
            imGradCost = exp( -imGaussGrad.^2 / (2 * sigmaGrad^2) );
            
            SecondOrderCliquePotential(:,2*i) = imGradCost(:);
            SecondOrderCliquePotential(:,2*i-1) = imGradCost(:);
            
            if flagDebugMode
                imseriesshow( imGradCost );
                set( gcf , 'Name' , sprintf( 'Second-Order Clique Potentials -- %s-Dir', directionLabels{i} ) );         
            end
            
        end
        
        SecondOrderCliquePotential = regularizationWeight * SecondOrderCliquePotential;
        
        
    %% Solve Maxflow using Kolmogorov-Boykov algorithm with a topology-preserving constaraint
    if flagDebugMode
        fprintf( 1 , '\n\nSolving the Maximum Flow Problem with Topology Constraints ... '  );
    end
    
    tic
    imSegMask = GcGridMaxflowTP( logical(imTopologyMask), FirstOrderCliquePotential_Fgnd, FirstOrderCliquePotential_Bgnd, SecondOrderCliquePotential, true );       
    imLabeledMask = bwlabeln(imSegMask, 2 * ndims(imInput));
    timeElapsed = toc;

    if flagDebugMode
        fprintf( 1 , 'It took %.2f seconds\n\n', timeElapsed  );
        
        imSegMaskRGB = label2rgbND(imLabeledMask);
        imseriesmaskshowrgb( imInput , imSegMaskRGB );

        title( 'Cell Segmentation Result using topology Preserving Graph-Cuts' , ...
               'FontSize' , 12 , 'FontWeight' , 'Bold' );    
        set( gcf , 'Name' , 'Cell Segmentation Result using topology Preserving Graph-Cuts: SegResult -- Topology Mask' );
        
    end
    
        
 end

