function [ imLabeledMask ] = segmentCellsUsingTopoGraphcut_ver2( imInput, varargin )
%% Segments cells in a uni/multi-cellular image using topology-preserving preserving graph cuts
%
%
    p = inputParser;
    ip.CaseSensitive = false;
    p.addRequired( 'imInput', @(x) (isnumeric(x) && ismember( ndims(x), [2,3] )) );
    p.parse( imInput );
    p.addParamValue( 'spacing', ones( 1, ndims(imInput) ), @(x) (isnumeric(x) && numel(x) == ndims(imInput)) );
    p.addParamValue( 'neighLabelRegularizationWeight', 1.0, @(x) (isscalar(x)) );    
    p.addParamValue( 'maxForegroundModes', 4, @(x) (isscalar(x)) );
    p.addParamValue( 'minCellDiameter', [], @(x) (isscalar(x)) ); % should be specified in physical space
    p.addParamValue( 'debugMode', true, @(x) (isscalar(x) && islogical(x)) );
    p.parse( imInput, varargin{:} );
    
    spacing = p.Results.spacing;
    neighLabelRegularizationWeight = p.Results.neighLabelRegularizationWeight;
    maxForegroundModes = p.Results.maxForegroundModes;
    minCellDiameter = p.Results.minCellDiameter;    
    flagDebugMode = p.Results.debugMode;  
    
    imsize = size(imInput);    
    imAdjusted = filterGaussND( imInput, 3.0, 'spacing', spacing );
    
    if flagDebugMode
        PrettryPrintStepDescription( sprintf( 'Segmenting cells using topology preserving graph-cuts in an image of size [ %s ]', ...
                                              sprintf(' %d ', size(imInput)) ));
    end
    
    %% Build Foreground and Background Pdf
    if flagDebugMode
        fprintf( '\n\nStep-1: Estimating foreground and background pdf ...\n\n' );
    end
    
    [imBackgroundMask, bgMean, bgSigma]  = estimateBackgroundArea( imAdjusted, ...
                                                                   'PostProcess', true );
    bgmodel = gmdistribution( bgMean, bgSigma^2 );    
    
    imForegroundMask = ~imBackgroundMask;
    fgmodel = cell(1,maxForegroundModes);
    AIC = zeros(1,4);
    numRandomFgPixels = min( 100000, 0.1 * numel(find(imForegroundMask)) );
    fgPixels = randsample( imAdjusted(imForegroundMask), numRandomFgPixels );
    for k = 1:maxForegroundModes
        fgmodel{k} = gmdistribution.fit( fgPixels, k, ...
                                         'Start', 'randSample', ...
                                         'Replicates', 3 );
        AIC(k) = fgmodel{k}.AIC;
    end
    [~,numComponents] = min(AIC);
    fgmodel = fgmodel{numComponents};
    
    if flagDebugMode
        
        fprintf( '\n%d foreground gaussians were found with mixture proportions - [ %s ]\n', ...
                 numComponents, ...
                 sprintf( ' %.3f ', fgmodel.PComponents ) );
             
        DisplayModelPdf( fgmodel , bgmodel );
        set( gcf, 'Name', 'Estimated foregrsound and background pdf' );
        
    end
        
    %% Generate Topology Mask   
    if flagDebugMode        
        fprintf( '\n\nStep-2: Estimating topology mask ...\n\n' );
    end
    
    if ~isempty( minCellDiameter )
        minCellDiameter = min( minCellDiameter ./ spacing ); % convert to image space
    end
    
    switch ndims(imInput)
       
        case 2            
            
            if isempty( minCellDiameter ) 
                
                % estimate object scale using connected component analysis
                imForegroundMask = imfill(imForegroundMask, 'holes');
                imDistMap = bwdist(~imForegroundMask);
                objScale = max( imDistMap(:) )/6.0;                
                
            else
                objScale = minCellDiameter/3.0;
            end            
            
            imTopologyMask = findCellTopology( imInput, ...
                                               'scale', objScale, 'plot', flagDebugMode );   
            imTopologyMask(~imForegroundMask) = 0;

        case 3        

            if isempty( minCellDiameter ) 
                
                % estimate object scale using connected component analysis
                % on the middle slice in the 2D + T matrix
                midSliceId = round( 0.5 * size(imInput,3) );
                imForegroundMask = imfill(imForegroundMask(:,:,midSliceId), 'holes');
                imDistMap = bwdist(~imForegroundMask);
                objScale = max( imDistMap(:) )/3.5;                
                
            else
                objScale = minCellDiameter/3.5;
            end                      
            
            % estimate the topology/numCells in 5 random frames and pick
            % the max            
            sliceList = randsample( 1:size(imInput,3), min(5,size(imInput,3)) ); 
            numComponents = zeros( size(sliceList) );
            for i = 1:numel(sliceList)
                mask = findCellTopology( imInput(:,:,sliceList(i)), 'scale', objScale );   
                mask(~imForegroundMask) = 0;
                L = bwlabel( mask > 0 );
                numComponents(i) = max( L(:) );
            end
            [numCells, id] = max( numComponents );
            maxComponentSlice = sliceList( id );
            
            mask = findCellTopology( imInput(:,:,maxComponentSlice), ...
                                     'scale', objScale, 'plot', flagDebugMode );            
            
            imTopologyMask = zeros( size(imInput) );            
            imTopologyMask(:,:,maxComponentSlice) = double( mask > 0 );
            
    end
    
    if flagDebugMode
        
        imseriesmaskshow( imInput, imTopologyMask );
        title( 'Topology Mask', 'FontSize' , 12 , 'FontWeight' , 'Bold' );
        set( gcf , 'Name' , 'Topology Mask' );        

    end
    

    %% Setup unary clique potentials of the MRF
    if flagDebugMode
        fprintf( '\n\nStep-3: Setting up unary clique potentials of the MRF energy ...\n\n' );
    end
    
        % Background
        FirstOrderCliquePotential_Bgnd = reshape( -log( pdf( bgmodel , imAdjusted(:) ) ) , imsize );

            % Fix Inf vals before normalizing
            maxBgndVal = max( FirstOrderCliquePotential_Bgnd( ~isinf( FirstOrderCliquePotential_Bgnd ) ) );
            FirstOrderCliquePotential_Bgnd( isinf( FirstOrderCliquePotential_Bgnd ) ) = maxBgndVal + 1;

            if flagDebugMode
                imseriesshow( FirstOrderCliquePotential_Bgnd );
                title( 'Unary Clique Potential -- Background', 'FontSize' , 12 , 'FontWeight' , 'Bold' );
                set( gcf , 'Name' , 'Unary Clique Potential -- Background' );         
            end
    
        % Foreground
        FirstOrderCliquePotential_Fgnd = reshape( -log( pdf( fgmodel , imAdjusted(:) ) ) , imsize );

            % Fix Inf vals before normalizing
            maxFgndVal = max( FirstOrderCliquePotential_Fgnd( ~isinf( FirstOrderCliquePotential_Fgnd ) ) );
            FirstOrderCliquePotential_Fgnd( isinf( FirstOrderCliquePotential_Fgnd ) ) = maxFgndVal + 1;

            if flagDebugMode
                imseriesshow( FirstOrderCliquePotential_Fgnd );
                title( 'Unary Clique Potential -- Foreground', 'FontSize' , 12 , 'FontWeight' , 'Bold' );
                set( gcf , 'Name' , 'Unary Clique Potential -- Foreground' );         
            end    
            
    %% Setup pairwise clique potentials of the MRF
    if flagDebugMode
        fprintf( '\n\nStep-4: Setting up pairwise-clique potentials of the MRF energy ...\n\n' );
    end
    
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
            sigmaGrad = std( imGaussGrad(imForegroundMask) );            
            imGradCost = exp( -imGaussGrad.^2 / (2 * sigmaGrad^2) );
            
            SecondOrderCliquePotential(:,2*i) = imGradCost(:);
            SecondOrderCliquePotential(:,2*i-1) = imGradCost(:);
            
            if flagDebugMode
                imseriesshow( imGradCost );
                set( gcf , 'Name' , sprintf( 'Second-Order Clique Potentials -- %s-Dir', directionLabels{i} ) );         
            end
            
        end
        
        SecondOrderCliquePotential = neighLabelRegularizationWeight * SecondOrderCliquePotential;
        
        
    %% Solve Maxflow using Kolmogorov-Boykov algorithm with a topology-preserving constaraint
    if flagDebugMode
        fprintf( 1 , '\n\nStep-5: Solving the Maximum Flow Problem with Topology Constraints ... '  );
    end
    
    tic
    imSegMask = GcGridMaxflowTP( logical(imTopologyMask), FirstOrderCliquePotential_Fgnd, FirstOrderCliquePotential_Bgnd, SecondOrderCliquePotential, true );       
    imLabeledMask = bwlabeln(imSegMask, 6);
    
    timeElapsed = toc;

    if flagDebugMode
        fprintf( 1 , 'It took %.2f seconds\n\n', timeElapsed  );
        
        imSegMaskRGB = label2rgbND(imLabeledMask);
        imseriesmaskshowrgb( imInput, imSegMaskRGB );

        title( 'Cell Segmentation Result using topology Preserving Graph-Cuts' , ...
               'FontSize' , 12 , 'FontWeight' , 'Bold' );    
        set( gcf , 'Name' , 'Cell Segmentation Result using topology Preserving Graph-Cuts: SegResult -- Topology Mask' );        
    end
    
        
 end

function PrettryPrintStepDescription( strStepDescription )

    strStar = strStepDescription;
    strStar(:) = '*';
    strStar = [ '****', strStar, '****' ];
    fprintf( '\n\n%s', strStar );
    fprintf( '\n    %s    ', strStepDescription );
    fprintf( '\n%s\n\n', strStar );
    
end

function DisplayModelPdf( fgmodel , bgmodel )

    % display foreground-background pdf
    if ~isempty( fgmodel ) && ~isempty( bgmodel )

        minMean = min( [ fgmodel.mu ; bgmodel.mu ] );
        maxMean = max( [ fgmodel.mu ; bgmodel.mu ] );
        minSigma = sqrt( min( [ fgmodel.Sigma(:) ; bgmodel.Sigma(:) ] ) );
        maxSigma = sqrt( max( [ fgmodel.Sigma(:) ; bgmodel.Sigma(:) ] ) );

        x = ( ( minMean - 3 * maxSigma ) : minSigma/20 : ( maxMean + 3 * maxSigma ) )';

        figure;

        hold on;

            plot( x , pdf( fgmodel , x ) , 'gx-' , 'LineWidth' , 2.0 );

            plot( x , pdf( bgmodel , x ) , 'rx-' , 'LineWidth' , 2.0 );

        hold off;

        title( 'Background and Foreground PDFs modeling using Mixture of Gaussians' , ...
               'FontSize' , 12 , 'FontWeight' , 'Bold' );

        set( gcf , 'Name' , 'Background and Foreground PDFs modeling using Mixture of Gaussians' );

        legend( { 'Foreground pdf' , 'Background pdf' } );    

    end    

end