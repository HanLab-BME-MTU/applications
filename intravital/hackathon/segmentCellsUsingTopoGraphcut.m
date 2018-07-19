function [ imLabeledMask ] = segmentCellsUsingTopoGraphcut( imInput, varargin )
%% Segments cells in a uni/multi-cellular image using topology-preserving preserving graph cuts
%
%
    p = inputParser;
    ip.CaseSensitive = false;
    p.addRequired( 'imInput', @(x) (isnumeric(x) && ismember( ndims(x), [2,3] )) );
    p.parse( imInput );
    p.addParamValue( 'spacing', ones( 1, ndims(imInput) ), @(x) (isnumeric(x) && numel(x) == ndims(imInput)) );
    p.addParamValue( 'neighLabelRegularizationWeight', 10.0, @(x) (isscalar(x)) );    
    p.addParamValue( 'maxForegroundModes', 3, @(x) (isscalar(x)) );
    p.addParamValue( 'debugMode', true, @(x) (isscalar(x) && islogical(x)) );
    p.parse( imInput, varargin{:} );
    
    spacing = p.Results.spacing;
    neighLabelRegularizationWeight = p.Results.neighLabelRegularizationWeight;
    maxForegroundModes = p.Results.maxForegroundModes;
    flagDebugMode = p.Results.debugMode;      
    
    imsize = size(imInput);    
    imInput = double( imInput );
    ImageIntensityRange = ComputeImageDynamicRange( imInput, 98.0 );        
    
    if flagDebugMode
        PrettryPrintStepDescription( sprintf( 'Segmenting cells using topology preserving graph-cuts in an image of size [ %s ]', ...
                                              sprintf(' %d ', size(imInput)) ));
    end    
    
    %% Prepare Topology Mask   
    if flagDebugMode        
        fprintf( '\n> Preparing topology mask ...\n\n' );
    end
    
    switch ndims(imInput)
       
        case 2
            
            [fgnd_seed_points, ~] = get_fgnd_bgnd_seeds_2d_points( imInput, ImageIntensityRange );

            imTopologyMask = zeros(size(imInput));    
            for i = 1:size(fgnd_seed_points,1)
                imTopologyMask( fgnd_seed_points(i,2), fgnd_seed_points(i,1) ) = i;
            end

        case 3        

            [fgnd_seed_points , ~] = get_fgnd_bgnd_seeds_2d_points( imInput(:,:,1), ImageIntensityRange );

            imTopologyMask = zeros( size(imInput) );
            for i = 1:size(fgnd_seed_points,1)
                imTopologyMask( fgnd_seed_points(i,2), fgnd_seed_points(i,1), 1) = i;
            end
            
    end
    
    if flagDebugMode
        fprintf( '\n\tThe image contains %d cells\n', size(fgnd_seed_points,1) );
    end                
    
    %% Setup and minimize MRF energy with a topology-preserving constaraint
    if flagDebugMode
        fprintf( '\n\n> Setting up and minimizing MRF energy with topology-preserving constraint...\n\n' );
        
        imBackgroundMaskForDisplay = zeros(size(imInput));
        imForegroundMaskForDisplay = zeros(size(imInput));        
        FirstOrderCliquePotential_Bgnd_ForDisplay = zeros(size(imInput));
        FirstOrderCliquePotential_Fgnd_ForDisplay = zeros(size(imInput));
        imEdgeCostForDisplay = zeros( [size(imInput), 2] );
        hModelPdf = figure;
        set( hModelPdf, 'Name', 'Estimated foreground and background pdf' );
    end
    
    imLabeledMask = zeros(size(imInput));
    
    tic
    for sliceId = 1:size(imInput,3)        
        
        if flagDebugMode && ndims(imInput) == 3
            fprintf( '\n\tSegmenting cells in slice %d/%d ...', sliceId, size(imInput,3) );
        end
        
        imCurSlice = mat2gray( imInput(:,:,sliceId) ) * 255.0; % standardize the intensity range 
        
        % do some denoising and smoothing
        imCurSlice = medfilt2( imCurSlice, [3,3] );
        imAdjusted = filterGaussND( imCurSlice, 1.5, 'spacing', spacing(1:2) );
        
        % prepare topology mask        
        if sliceId > 1
           
            topologyPropagationType = 'prevSegMaskCentroid';
            switch topologyPropagationType

                case 'prevSegMaskShrinked'
                    
                    % set topology as a shrinked version of the segmentation in previous slice     
                    imCurTopologyMask = zeros( size(imCurSlice) );
                    imPrevLabeledSedMask = imLabeledMask(:,:,sliceId-1);
                    for cid = 1:size(fgnd_seed_points,1)          
                        imShrinkedMask = bwmorph( imPrevLabeledSedMask == cid, 'shrink', 3 );
                        imCurTopologyMask( imShrinkedMask ) = cid;
                    end   

                    imTopologyMask(:, :, sliceId) = imCurTopologyMask;
                    
                case 'prevSegMaskCentroid'
                    
                    imCurTopologyMask = zeros( size(imCurSlice) );
                    imPrevLabeledSedMask = imLabeledMask(:,:,sliceId-1);
                    for cid = 1:size(fgnd_seed_points,1)          
                        [yind, xind] = ind2sub( size(imCurTopologyMask), find(imPrevLabeledSedMask == cid) );
                        ptCentroid = round( mean( [xind, yind], 1 ) );
                        imCurTopologyMask( ptCentroid(2), ptCentroid(1) ) = cid;
                    end   

                    imTopologyMask(:, :, sliceId) = imCurTopologyMask;
                    
            end
            
        else            
            imCurTopologyMask = imTopologyMask(:, :, sliceId);            
        end
        
        
        % Estimate foreground and background pdf
        if sliceId == 1
            
            [imBackgroundMask, bgMean, bgSigma]  = estimateBackgroundArea( imAdjusted, ...
                                                                           'PostProcess', true );
            bgmodel = gmdistribution( bgMean, bgSigma^2 );    


            imForegroundMask = imAdjusted > thresholdFluorescenceImage( imAdjusted );
            fgmodel = cell(1,maxForegroundModes);
            BIC = zeros(1,maxForegroundModes);
            numRandomFgPixels = round( min( 100000, 0.1 * numel(find(imForegroundMask)) ) );
            fgPixels = randsample( imAdjusted(imForegroundMask), numRandomFgPixels );
            for k = 1:maxForegroundModes
                fgmodel{k} = gmdistribution.fit( fgPixels, k, ...
                                                 'Start', 'randSample', ...
                                                 'Replicates', 3 );
                BIC(k) = fgmodel{k}.BIC;
            end
            [~,numForegroundModes] = min(BIC);
            fgmodel = fgmodel{numForegroundModes};

        else
            
            [imBackgroundMask, bgMean, bgSigma]  = estimateBackgroundArea( imAdjusted, ...
                                                                           'PostProcess', true );
            bgmodel = gmdistribution( bgMean, bgSigma^2 );    


            imForegroundMask = imAdjusted > thresholdFluorescenceImage( imAdjusted );
            
            numRandomFgPixels = round( min( 100000, 0.1 * numel(find(imForegroundMask)) ) );
            fgPixels = randsample( imAdjusted(imForegroundMask), numRandomFgPixels );
            initModelParams = struct;
            initModelParams.mu = fgmodel.mu;
            initModelParams.Sigma = fgmodel.Sigma;
            initModelParams.PComponents = fgmodel.PComponents;
            fgmodel = gmdistribution.fit( fgPixels, numForegroundModes, ...
                                          'Start', initModelParams );
            
        end       

        if flagDebugMode
            
            DisplayModelPdf( fgmodel , bgmodel, hModelPdf );
            
            fprintf( '\n\n\tforeground was modelled with %d-gaussian mixture model:\n\t mean - [ %s ]\n\t stddev - [ %s ] \n\t proportions - [ %s ]\n\n', ...
                     numForegroundModes, ...
                     sprintf( ' %.3f ', fgmodel.mu ), ...
                     sprintf( ' %.3f ', sqrt(fgmodel.Sigma) ), ...
                     sprintf( ' %.3f ', fgmodel.PComponents ) );
            
            imBackgroundMaskForDisplay(:,:,sliceId) = imBackgroundMask;
            imForegroundMaskForDisplay(:,:,sliceId) = imForegroundMask;
        end
        
        % Setup unary clique potentials of the MRF

            % Background
            FirstOrderCliquePotential_Bgnd = reshape( -log( pdf( bgmodel , imAdjusted(:) ) ) , imsize(1:2) );

                % Fix Inf vals before normalizing
                maxBgndVal = max( FirstOrderCliquePotential_Bgnd( ~isinf( FirstOrderCliquePotential_Bgnd ) ) );
                FirstOrderCliquePotential_Bgnd( isinf( FirstOrderCliquePotential_Bgnd ) ) = maxBgndVal + 1;
                FirstOrderCliquePotential_Bgnd = FirstOrderCliquePotential_Bgnd - min( FirstOrderCliquePotential_Bgnd(:) );
                
                if flagDebugMode
                    FirstOrderCliquePotential_Bgnd_ForDisplay(:,:,sliceId) = FirstOrderCliquePotential_Bgnd;
                end

            % Foreground
            FirstOrderCliquePotential_Fgnd = reshape( -log( pdf( fgmodel , imAdjusted(:) ) ) , imsize(1:2) );

                % Fix Inf vals before normalizing
                maxFgndVal = max( FirstOrderCliquePotential_Fgnd( ~isinf( FirstOrderCliquePotential_Fgnd ) ) );
                FirstOrderCliquePotential_Fgnd( isinf( FirstOrderCliquePotential_Fgnd ) ) = maxFgndVal + 1;
                FirstOrderCliquePotential_Fgnd = FirstOrderCliquePotential_Fgnd - min( FirstOrderCliquePotential_Fgnd(:) );
                
                if flagDebugMode
                    FirstOrderCliquePotential_Fgnd_ForDisplay(:,:,sliceId) = FirstOrderCliquePotential_Fgnd;
                end       
                
        % Setup pairwise clique potentials of the MRF
        SecondOrderCliquePotential = zeros( numel(imCurSlice), 2*ndims(imCurSlice) );
        pairwisePotentialtype = 'steerableFilterResponse';
        
        switch pairwisePotentialtype
            
            case 'gradient'
                
                for i = 1:2

                    % edge cost
                    imGaussGrad = filterGaussGradND(imCurSlice, 1.5, i, 'spacing', spacing(1:2) );           
                    sigmaGrad = std( imGaussGrad(imBackgroundMask) ); 
                    meanGrad = mean( imGaussGrad(imBackgroundMask) ); 
                    imEdgeCost = exp( -(imGaussGrad - meanGrad).^2 / (2 * sigmaGrad^2) );

                    SecondOrderCliquePotential(:,2*i) = imEdgeCost(:) + 1;
                    SecondOrderCliquePotential(:,2*i-1) = imEdgeCost(:) + 1;

                    if flagDebugMode
                        imEdgeCostForDisplay(:,:,sliceId,i) = imEdgeCost;
                    end

                end
                
            case 'steerableFilterNms'

                [res, theta, nms,~] = steerableDetector(imCurSlice, 3, 2);                
                
                orientationVec = cell(1,2);
                orientationVec{1} = cos( pi/2.0 + theta );
                orientationVec{2} = sin( pi/2.0 + theta );
                
                for i = 1:2

                    % edge cost   
                    curNms = (nms .* orientationVec{i}).^2;
                    sigmaNms = std( curNms(imBackgroundMask) ); 
                    meanNms = mean( curNms(imBackgroundMask) ); 
                    imEdgeCost = exp( -(curNms - meanNms).^2 / (2 * sigmaNms^2) );

                    SecondOrderCliquePotential(:,2*i) = imEdgeCost(:) + 1;
                    SecondOrderCliquePotential(:,2*i-1) = imEdgeCost(:) + 1;

                    if flagDebugMode
                        imEdgeCostForDisplay(:,:,sliceId,i) = imEdgeCost;
                    end

                end

            case 'steerableFilterResponse'

                [res, theta, nms,~] = steerableDetector(imCurSlice, 3, 2);                
                
                orientationVec = cell(1,2);
                orientationVec{1} = cos( pi/2.0 + theta );
                orientationVec{2} = sin( pi/2.0 + theta );
                
                for i = 1:2

                    % edge cost   
                    curRes = (res .* orientationVec{i}).^2;
                    sigmaRes = std( curRes(imBackgroundMask) ); 
                    meanRes = mean( curRes(imBackgroundMask) ); 
                    imEdgeCost = exp( -(curRes - meanRes).^2 / (2 * sigmaRes^2) );

                    SecondOrderCliquePotential(:,2*i) = imEdgeCost(:) + 1;
                    SecondOrderCliquePotential(:,2*i-1) = imEdgeCost(:) + 1;

                    if flagDebugMode
                        imEdgeCostForDisplay(:,:,sliceId,i) = imEdgeCost;
                    end

                end
                
        end                

        SecondOrderCliquePotential = neighLabelRegularizationWeight * SecondOrderCliquePotential;

        % Solve Maxflow using Kolmogorov-Boykov algorithm with a topology-preserving constaraint        
        imCurSegMask = GcGridMaxflowTP( logical(imCurTopologyMask), ...
                                        FirstOrderCliquePotential_Fgnd, ...
                                        FirstOrderCliquePotential_Bgnd, ...
                                        SecondOrderCliquePotential, true );       
        
        L = bwlabeln(imCurSegMask, 4);
        imCurLabeledMask = zeros( size(imCurSlice) );        
        stats = regionprops( imCurTopologyMask, 'Centroid' );
        for cid = 1:size(fgnd_seed_points,1) 
            curCentoid = round( stats(cid).Centroid );
            seedLabel = L( curCentoid(2), curCentoid(1) );
            imCurLabeledMask( L == seedLabel ) = cid;
        end                
        imLabeledMask(:,:,sliceId) = imCurLabeledMask;
        
    end
    timeElapsed = toc;   

    if flagDebugMode
        
        fprintf( 1 , '\n\tThe segmentation took %.2f seconds\n\n', timeElapsed  );
        
        % display topology mask
        if strcmp(topologyPropagationType, 'prevMaskCentroid')
            imTopologyMask = imdilate( imTopologyMask, ones(10,10) );
        else
            imTopologyMask(:,:,1) = imdilate( imTopologyMask(:,:,1), ones(10,10) );
        end
        
        imTopologyMaskRGB = label2rgbND( imTopologyMask );
        imseriesmaskshowrgb( imInput, imTopologyMaskRGB );
        title( 'Topology Mask', 'FontSize' , 12 , 'FontWeight' , 'Bold' );
        set( gcf , 'Name' , 'Topology Mask' );        
        
        % display foreground-background masks used to build appearance models
        imseriesmaskshow( imInput, {imForegroundMaskForDisplay, imBackgroundMaskForDisplay} );
        set( gcf , 'Name' , 'Foreground-background masks used for appearance models' );         
        
        % display background unary potential
        imseriesshow( FirstOrderCliquePotential_Bgnd_ForDisplay );
        title( 'Unary Clique Potential -- Background', 'FontSize' , 12 , 'FontWeight' , 'Bold' );
        set( gcf , 'Name' , 'Unary Clique Potential -- Background' );         
        
        % display foreground unary potential
        imseriesshow( FirstOrderCliquePotential_Fgnd_ForDisplay );
        title( 'Unary Clique Potential -- Foreground', 'FontSize' , 12 , 'FontWeight' , 'Bold' );
        set( gcf , 'Name' , 'Unary Clique Potential -- Foreground' );         
        
        % display pairwise potentisls
        directionLabels = { 'Y', 'X' };
        
        for i = 1:2
            imseriesshow( imEdgeCostForDisplay(:,:,:,i) );
            set( gcf , 'Name' , sprintf( 'Second-Order Clique Potentials -- %s-Dir', directionLabels{i} ) );         
        end
        
        % display segmentation result
        imSegMaskRGB = label2rgbND(imLabeledMask);
        imseriesmaskshowrgb( imInput, imSegMaskRGB );
        title( 'Cell Segmentation Result using topology Preserving Graph-Cuts' , ...
               'FontSize' , 12 , 'FontWeight' , 'Bold' );    
        set( gcf , 'Name' , 'Cell Segmentation Result using topology Preserving Graph-Cuts: SegResult' );        
        
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

function DisplayModelPdf( fgmodel , bgmodel, figHandle )

    % display foreground-background pdf
    if ~isempty( fgmodel ) && ~isempty( bgmodel )

        minMean = min( [ fgmodel.mu ; bgmodel.mu ] );
        maxMean = max( [ fgmodel.mu ; bgmodel.mu ] );
        minSigma = sqrt( min( [ fgmodel.Sigma(:) ; bgmodel.Sigma(:) ] ) );
        maxSigma = sqrt( max( [ fgmodel.Sigma(:) ; bgmodel.Sigma(:) ] ) );

        x = ( ( minMean - 3 * maxSigma ) : minSigma/20 : ( maxMean + 3 * maxSigma ) )';

        if ~exist( 'figHandle', 'var' )
            figure;
        end

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