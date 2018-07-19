clc
clear 
close all

datadim = 3;

overlap = 0.2;
theta = (0:120:359);
fgMeanVar = [ 200, 20 ];
bgMeanVar = [ 180, 20 ];

switch datadim
    
    case 2
    
        im = zeros(200,200);
        cellRadius = [10, 20];

        imsize = size(im);
        minObjectDiameter = 2 * min(cellRadius);

        fgMeanVar = [ 200, 20 ];
        bgMeanVar = [ 180, 20 ];

        fgGmObj = gmdistribution( fgMeanVar(:,1), reshape(fgMeanVar(:,2), [1,1,size(fgMeanVar,1)]) );
        bgGmObj = gmdistribution( bgMeanVar(:,1), reshape(bgMeanVar(:,2), [1,1,size(bgMeanVar,1)]) );

        [X,Y] = meshgrid(1:imsize(1),1:imsize(2));

        im(:) = random(bgGmObj,numel(im));
        fgMask = zeros( size(im) );
        %im(:) = bgMeanVar(1);

        for i = 1:numel(theta)    
            xc = 0.5 * imsize(2) + max(cellRadius) * (1 - overlap) * cosd( theta(i) );
            yc = 0.5 * imsize(1) + max(cellRadius) * (1 - overlap) * sind( theta(i) );

            pts = [ X(:) -  xc, Y(:) - yc ];
            pts = pts * [ cosd(theta(i)) -sind(theta(i)); sind(theta(i)), cosd(theta(i)) ];

            curEllipseInd = find( (pts(:,1).^2 / (cellRadius(2))^2)  + (pts(:,2).^2 / (cellRadius(1))^2) - 1 <= 0 );    

            % randomly select a foreground distrubution and draw samples from it
            im( curEllipseInd ) = random(fgGmObj,numel(curEllipseInd)); 
            %im( curEllipseInd ) = 1; 
            fgMask( curEllipseInd ) = 1;

        end        
        
    case 3
        
        im = zeros(100,100,15);
        cellRadius = [5, 10, 5];

        imsize = size(im);
        minObjectDiameter = 2 * min(cellRadius);

        fgGmObj = gmdistribution( fgMeanVar(:,1), reshape(fgMeanVar(:,2), [1,1,size(fgMeanVar,1)]) );
        bgGmObj = gmdistribution( bgMeanVar(:,1), reshape(bgMeanVar(:,2), [1,1,size(bgMeanVar,1)]) );

        [X,Y,Z] = meshgrid(1:imsize(1),1:imsize(2),1:imsize(3));

        im(:) = random(bgGmObj,numel(im));
        fgMask = zeros( size(im) );
        %im(:) = bgMeanVar(1);

        for i = 1:numel(theta)    
            xc = 0.5 * imsize(2) + max(cellRadius) * (1 - overlap) * cosd( theta(i) );
            yc = 0.5 * imsize(1) + max(cellRadius) * (1 - overlap) * sind( theta(i) );
            zc = 0.5 * imsize(3);

            pts = [ X(:) -  xc, Y(:) - yc ];
            pts = pts * [ cosd(theta(i)) -sind(theta(i)); sind(theta(i)), cosd(theta(i)) ];
            pts(:,3) = Z(:) - zc;

            curEllipseInd = find( (pts(:,1).^2 / (cellRadius(2))^2)  + (pts(:,2).^2 / (cellRadius(1))^2) + (pts(:,3).^2 / (cellRadius(3))^2) - 1 <= 0 );    

            % randomly select a foreground distrubution and draw samples from it
            im( curEllipseInd ) = random(fgGmObj,numel(curEllipseInd)); 
            %im( curEllipseInd ) = 1; 
            fgMask( curEllipseInd ) = 1;

        end
        
end

% compare different methods for cells seed point detection   
meanCellDiameter = 2 * mean( cellRadius );    
cellDiameterRange = 2 * ( mean(cellRadius) + [-2,2] );
numScales = 5;
imCellSeedDetectionResult = {};
strSeedDetectionAlgorithm = {};

    % Local Intensity Maxima in Blurred Image
    fprintf( '\n\nDetecting seed points as local maxima in Blurred Image ...\n\n' );

    strSeedDetectionAlgorithm{end+1} = 'IntensityMaxima';
    imCellSeedDetectionResult{end+1} = detect_cell_seeds_IntensityMaxima( im, ...
                                                                          meanCellDiameter, ...
                                                                          'debugMode', true );

    % Local Maxima in LoG Filtered Image
    fprintf( '\n\nDetecting seed points as local maxima in LoG Filtered Image ...\n\n' );

    strSeedDetectionAlgorithm{end+1} = 'FixedScaleLoG';
    imCellSeedDetectionResult{end+1} = detectBlobsUsingLoG( im, ...
                                                            meanCellDiameter, ...
                                                            'debugMode', true );

    % Local Maxima in a Multiscale LoG Filtered Image
    fprintf( '\n\nDetecting seed points as local maxima in Multiscale LoG Filtered Image ...\n\n' );

    strSeedDetectionAlgorithm{end+1} = 'MultiScaleLoG';
    
    imCellSeedDetectionResult{end+1} = detectBlobsUsingMultiscaleLoG( im, ...
                                                                      cellDiameterRange, ...
                                                                      'numLoGScales', numScales, ...
                                                                      'debugMode', true );

    % Local Maxima in a Multiscale LoG Filtered Image
    fprintf( '\n\nDetecting seed points as local maxima in Adaptive Multiscale LoG Filtered Image ...\n\n' );

    strSeedDetectionAlgorithm{end+1} = 'AdaptiveMultiScaleLoG';
    imCellSeedDetectionResult{end+1} = detect_cell_seeds_adaptive_multiscale_LoG( im, ...
                                                                                  fgMask, ...                                                                                  
                                                                                  'cellDiameterRange', cellDiameterRange, ...
                                                                                  'numLoGScales', numScales, ...
                                                                                  'debugMode', true );

    % Local Maxima in a Multiscale LoBG Filtered Image
    fprintf( '\n\nDetecting seed points as local maxima in Adaptive Multiscale LoG Filtered Image ...\n\n' );

    strSeedDetectionAlgorithm{end+1} = 'MultiScaleLoBG';
    imCellSeedDetectionResult{end+1} = detectBlobsUsingMultiscaleLoBG( im, ...
                                                                       cellDiameterRange, ...
                                                                       'numLoGScales', numScales, ...
                                                                       'debugMode', true );

    % Local Maxima in a radial symmetry transform
    fprintf( '\n\nDetecting seed points as local maxima in multiscale radial symmetry transform ...\n\n' );

    strSeedDetectionAlgorithm{end+1} = 'RadialSymmetryTransform';
    roiMask = imdilate( fgMask, ones(3,3) );
    imCellSeedDetectionResult{end+1} = detect_cell_seeds_radial_symmetry( im, ...
                                                                          0.5 * cellDiameterRange, ...
                                                                          'roiMask', roiMask, ...
                                                                          'numRadiusSamples', numScales, ...
                                                                          'flagParallelize', false, ...
                                                                          'debugMode', true );

    % suppress seed points outside the thresholded foreground region
    for i = 1:numel( imCellSeedDetectionResult )
        imCellSeedDetectionResult{i}( ~fgMask ) = 0;
    end

    % display result
    imSeedMaskForDisplay = imCellSeedDetectionResult;
    for i = 1:numel( imSeedMaskForDisplay )
        imSeedMaskForDisplay{i} = imdilate( imSeedMaskForDisplay{i}, strel('disk', 3) );
    end
    imseriesmaskshow( im, imSeedMaskForDisplay );
    set( gcf, 'Name', [ 'Cell Seed Point Detection Results', ...
                        sprintf( ' -- %s', strSeedDetectionAlgorithm{:} ) ] );
