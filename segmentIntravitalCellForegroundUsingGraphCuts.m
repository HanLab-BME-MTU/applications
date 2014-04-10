function [ imCellForegroundMask ] = segmentIntravitalCellForegroundUsingGraphCuts( im, varargin )

    p = inputParser;
    p.addRequired( 'im', @(x) (isnumeric(x) && ismember( ndims(x), [2,3] )) );
    p.parse( im );    
    p.addParamValue( 'spacing', ones( 1, ndims(im) ), @(x) (isnumeric(x) && numel(x) == ndims(im)) );    
    p.addParamValue( 'regularizationWeight', 1.0, @(x) (isscalar(x)) );
    p.addParamValue( 'downsamplingFactor', 1, @(x) (isscalar(x) || (isnumeric(x) && numel(x) == ndims(im))) );
    p.addParamValue( 'debugMode', true, @(x) (isscalar(x) && islogical(x)) );    
    p.addParamValue( 'flagParallelize', false, @(x) (isscalar(x) && islogical(x)) );
    p.parse( im, varargin{:} );

    spacing = p.Results.spacing;    
    regularizationWeight = p.Results.regularizationWeight;
    downsamplingFactor = p.Results.downsamplingFactor;    
    
    flagDebugMode = p.Results.debugMode;    
    flagParallelize = p.Results.flagParallelize;
    
    % downsample the image if requested (this speeds up things)
    imUnsampled = im;
    if ~(isscalar( downsamplingFactor ) && downsamplingFactor == 1)
        flagDownsampled = true;
        im = imresizend( im, downsamplingFactor );      
        spacing = spacing ./ downsamplingFactor;
    end        
        
    imsize = size( im );

    % open matlab pool if parallelization is requested
    if flagParallelize
        flagPoolOpenedAlready = matlabpool( 'size' ) > 0;        
        if ~flagPoolOpenedAlready 
            matlabpool open;
        end            
    end
    
    % pre-processing
    if flagDebugMode
        fprintf( '\n\nPre-processing the image ... \n\n'  );
    end
    
    ImageIntensityRange = ComputeImageDynamicRange( im, 98.0 );
    imAdjusted = mat2gray( im, ImageIntensityRange ) * 4096;    
    imAdjusted = matitk( 'FMEDIAN', [2,2,2], imAdjusted );    
    
    % Compute high-confident foreground and background regions
    if flagDebugMode
        fprintf( '\n\nComputing high-confidence foreground and background regions ... \n\n'  );
    end
    
        % log transform        
        imAdjustedLog = mat2gray( log( 1 + imAdjusted ) );
        threshGlobalVal = thresholdOtsu( imAdjustedLog );      
        
        % two-levels of otsu thresholding in each slice
        localThresholdFactor = 0.70;
        imThresh = zeros( size( imAdjusted ) );
        imThreshHighConfidence = zeros( size( imAdjusted ) );

        if flagParallelize
            parfor sliceId = 1:imsize(3)     

                imSlice = imAdjustedLog(:,:,sliceId);                                                                 
                curThreshVal = thresholdOtsu( imSlice );                    
                if curThreshVal < (localThresholdFactor * threshGlobalVal)
                    continue;
                else
                    imThresh(:,:,sliceId) = (imSlice > curThreshVal);
                    curThreshValHighConfidence = thresholdOtsu( imSlice(imSlice > curThreshVal) );
                    imThreshHighConfidence(:,:,sliceId) = (imSlice > curThreshValHighConfidence);
                end                    

            end 
        else        
            for sliceId = 1:imsize(3)      

                imSlice = imAdjustedLog(:,:,sliceId);                                                                 
                curThreshVal = thresholdOtsu(imSlice);                    
                if curThreshVal < (localThresholdFactor * threshGlobalVal)
                    continue;
                else
                    imThresh(:,:,sliceId) = (imSlice > curThreshVal);
                    curThreshValHighConfidence = thresholdOtsu( imSlice(imSlice > curThreshVal) );
                    imThreshHighConfidence(:,:,sliceId) = (imSlice > curThreshValHighConfidence);
                end                    

            end                
        end
        
        % compute high-confidence foreground and background seed regions
        imForegroundSeedRegion = imThreshHighConfidence;
        imBackgroundSeedRegion = imerode( 1 - imThresh, streldisknd([2,2,1]) );
        
        % Display the computed seed regions
        if flagDebugMode        
            imseriesmaskshow( im, { imForegroundSeedRegion, imBackgroundSeedRegion } );
            set( gcf, 'Name', 'Highly Confident Foreground and Background Seed Regions' );
        end
    
    % Setup unary clique potentials
    infval = 100000;
    
        % Background
        FirstOrderCliquePotential_Bgnd = zeros( imsize );
        FirstOrderCliquePotential_Bgnd( imForegroundSeedRegion > 0 ) = infval;
        
        % Foreground
        FirstOrderCliquePotential_Fgnd = zeros( imsize );
        FirstOrderCliquePotential_Fgnd( imBackgroundSeedRegion > 0 ) = infval;
        
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
        
        imGrad = cell(1,ndims(im));
        spacingVec = num2cell( spacing );
        [imGrad{:}] = gradient( imAdjusted, spacingVec{:} );
        imGrad = imGrad( [2,1,3:ndims(im)] );
        
        imHomogenousRegion = (imBackgroundSeedRegion > 0 | imForegroundSeedRegion > 0);
        
        for i = 1:ndims(im)
            
            % edge cost
            %imGaussGrad = filterGaussGradND(imAdjusted, min(spacing), i, 'spacing', spacing );
            %sigmaGrad = std( imGaussGrad(imBackgroundSeedRegion > 0) );            
            sigmaGrad = std( imGrad{i}(imHomogenousRegion > 0) );   
            meanGrad = mean( imGrad{i}(imHomogenousRegion > 0) );   
            imGradCost = exp( -(imGrad{i} - meanGrad).^2 / (2 * sigmaGrad^2) );
            imGradCost( imGrad{i} < meanGrad ) = 1;
            
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
                imseriesshow( imGrad{i} );
                set( gcf , 'Name' , sprintf( 'Gradient -- %s-Dir', directionLabels{i} ) );         
                
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
        fprintf( 1 , '\n\nSegmenting cell foreground using MRF ... '  );
    end

        tic
        [ energy , seg_labels ] = MinimizeBinaryOrder2MRFEnergy_GraphCuts( FirstOrderCliquePotential, SecondOrderCliques, SecondOrderCliquePotential );       
        timeElapsed = toc;

        if flagDebugMode
            fprintf( 1 , 'It took %.2f seconds\n\n', timeElapsed  );
        end
        
        imCellForegroundMask = reshape( seg_labels , imsize );
        imCellForegroundMask = double( imCellForegroundMask > 0 );

   % Display the result
   if flagDebugMode          

       imseriesmaskshow( im, { imCellForegroundMask, imForegroundSeedRegion, imBackgroundSeedRegion } );

       title( 'Seg Result overlayed with foreground-background seed region' , ...
              'FontSize' , 12 , 'FontWeight' , 'Bold' );    
       set( gcf , 'Name' , 'SegResult -- foreground seed region --- background seed region' );

   end        
        
    % upsampled result if it was downsampled
    if flagDownsampled
        imCellForegroundMask = imresizend( imCellForegroundMask, size(imUnsampled) ./ size(im), ...
                                           'interpolationMethod', 'nearest' );                                             
    end
    
    % post-process it to remove any crap and smooth the boundaries a bit
    imCellForegroundMask = imopen( imCellForegroundMask, streldisknd([3,3,1]) );
    imCellForegroundMask = imfill( imCellForegroundMask );        
        
   % Display the result
   if flagDebugMode           

       imseriesmaskshow( imUnsampled, { imCellForegroundMask } );

       title( 'Cell Foreground Segmentation using MRF' , ...
              'FontSize' , 12 , 'FontWeight' , 'Bold' );    
       set( gcf , 'Name' , 'Cell Foreground Segmentation using MRF' );

   end        
    
end