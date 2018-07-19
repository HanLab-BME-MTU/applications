
clc
clear
close all

% generate image
image_type = 'square_noisy_3d';

switch image_type 
    
    case 'flat_square_2d'
        
        imInput = zeros( 512 , 512 );
        image_size = size( imInput );
        imInput( image_size(1) * 0.25 : image_size(1) * 0.75 , image_size(2) * 0.25 : image_size(2) * 0.75 ) = 255;

    case 'square_weakbnd_2d'        
        
        imInput = zeros( 512 , 512 );
        image_size = size( imInput );
        dvec = ones( image_size(1), 1 );
        gap = 0.05;
        dvec( round(0.5 * image_size(1) + (-gap*image_size(1):gap*image_size(1))) ) = 0;        
        imInput( diag(dvec) > 0 ) = 1;

        fgnd_seed_points{1} = round( 0.5 * image_size + 0.4 * image_size .* [-1, 1] );
        fgnd_seed_points{2} = round( 0.5 * image_size + 0.4 * image_size .* [1, -1] );
        
    case 'square_weakbnd_random_2d'        
        
        imInput = zeros( 512 , 512 );
        image_size = size( imInput );
        dvec = ones( image_size(1), 1 );
        rperm = randperm( numel(dvec) );
        fracDrop = 0.5;
        dvec( rperm(1:round(fracDrop*numel(dvec))) ) = 0;
        imInput = diag(dvec) * 255;
        
        fgnd_seed_points{1} = round( 0.5 * image_size + 0.4 * image_size .* [-1, 1] );
        fgnd_seed_points{2} = round( 0.5 * image_size + 0.4 * image_size .* [1, -1] );
        
    case 'square_noisy_2d'

        imInput = zeros( 512 , 512 );

        gmobj_bgnd_ex = gmdistribution( 50 , 30 );
        imInput_bgnd_ex = reshape( random( gmobj_bgnd_ex , numel( imInput ) ) , size( imInput ) );

        gmobj_fgnd_ex = gmdistribution( 90 , 30 );
        imInput_fgnd_ex = reshape( random( gmobj_fgnd_ex , numel( imInput ) ) , size( imInput ) );

        imInput = imInput_bgnd_ex;

        image_size = size( imInput );

        imInput( image_size(1) * 0.25 : image_size(1) * 0.75 , image_size(2) * 0.25 : image_size(2) * 0.75 ) = imInput_fgnd_ex( image_size(1) * 0.25 : image_size(1) * 0.75 , image_size(2) * 0.25 : image_size(2) * 0.75 );

        imInput = round( imInput );

        
    case 'ellipse_noisy_2d'

        imInput = zeros( 512 , 512 );

        gmobj_bgnd_ex = gmdistribution( 50 , 30 );
        imInput_bgnd_ex = reshape( random( gmobj_bgnd_ex , numel( imInput ) ) , size( imInput ) );

        gmobj_fgnd_ex = gmdistribution( 90 , 30 );
        imInput_fgnd_ex = reshape( random( gmobj_fgnd_ex , numel( imInput ) ) , size( imInput ) );

        imInput = imInput_bgnd_ex;

        image_size = size( imInput );
        majorRad = 0.35 * image_size(2);
        minorRad = 0.10 * image_size(1);
        yc = round(0.5 * image_size(1));
        xc = round(0.5 * image_size(2));
        [X,Y] = meshgrid(1:image_size(2),1:image_size(1));
        pts = [ X(:) - xc, Y(:) - yc ];
        
        ellipsePixInd = find( (pts(:,1).^2 / majorRad^2)  + (pts(:,2).^2 / minorRad^2) - 1 <= 0 );
        
        imInput(ellipsePixInd) = imInput_fgnd_ex(ellipsePixInd);    
        
        gap = 0.01;
        imInput( round(yc + (-minorRad : minorRad)), xc ) = 0;
        imInput( round(yc + (-gap * 2 * minorRad : gap * 2 * minorRad)), xc ) = imInput_fgnd_ex( round(yc + (-gap * 2 * minorRad : gap * 2 * minorRad)), xc );
        
    case 'ellipse_noisy_2d_interactive'

        imInput = zeros( 512 , 512 );

        gmobj_bgnd_ex = gmdistribution( 50 , 30 );
        imInput_bgnd_ex = reshape( random( gmobj_bgnd_ex , numel( imInput ) ) , size( imInput ) );

        gmobj_fgnd_ex = gmdistribution( 90 , 30 );
        imInput_fgnd_ex = reshape( random( gmobj_fgnd_ex , numel( imInput ) ) , size( imInput ) );

        imInput = imInput_bgnd_ex;

        image_size = size( imInput );
        majorRad = 0.35 * image_size(2);
        minorRad = 0.10 * image_size(1);
        yc = round(0.5 * image_size(1));
        xc = round(0.5 * image_size(2));
        [X,Y] = meshgrid(1:image_size(2),1:image_size(1));
        pts = [ X(:) - xc, Y(:) - yc ];
        
        ellipsePixInd = find( (pts(:,1).^2 / majorRad^2)  + (pts(:,2).^2 / minorRad^2) - 1 <= 0 );
        
        imInput(ellipsePixInd) = imInput_fgnd_ex(ellipsePixInd);
        
       
    case 'square_noisy_3d'

        imInput = zeros( 100, 100 );

        gmobj_bgnd_ex = gmdistribution( 50 , 30 );
        imInput_bgnd_ex = reshape( random( gmobj_bgnd_ex , numel( imInput ) ) , size( imInput ) );

        gmobj_fgnd_ex = gmdistribution( 90 , 30 );
        imInput_fgnd_ex = reshape( random( gmobj_fgnd_ex , numel( imInput ) ) , size( imInput ) );

        imInput = imInput_bgnd_ex;

        image_size = size( imInput );

        imInput( round(image_size(1) * 0.25 : image_size(1) * 0.75) , ...
                 round(image_size(2) * 0.25 : image_size(2) * 0.75) ) = imInput_fgnd_ex( round(image_size(1) * 0.25 : image_size(1) * 0.75), round(image_size(2) * 0.25 : image_size(2) * 0.75) );

        imInput = round( imInput );      
        imInput = repmat(imInput, [1, 1, 10]);
end

% get seed points
if ~exist( 'fgnd_seed_points', 'var' )
    
    switch ndims(imInput)

        case 2

            [ fgnd_seed_points , bgnd_seed_points ] = get_fgnd_bgnd_seeds_2d_strokes( imInput );

        case 3

            [ fgnd_seed_points , bgnd_seed_points ] = get_fgnd_bgnd_seeds_3d_strokes( imInput );

    end
    
end

SEED_NEIGH = 3;
objSeedMask = zeros(size(imInput));
for i = 1:numel(fgnd_seed_points)
    curSeedIndices = sub2ind(size(imInput), fgnd_seed_points{i}(:,2), fgnd_seed_points{i}(:,1)); 
    objSeedMask(curSeedIndices) = i;
end
objSeedMask = imdilate( objSeedMask, streldisknd(SEED_NEIGH*ones(1,ndims(imInput))) );
objSeedIndices = find(objSeedMask);

seedIndices = objSeedIndices;
seedLabels = objSeedMask(objSeedIndices) - 1;
numLabels = numel(fgnd_seed_points);

% run random walker segmentation algorithm
tic

[segMask, pixToLabelProbabilityMap] = random_walker_segmentation( imInput, seedIndices, seedLabels );
%[segMask, pixToLabelProbabilityMap] = random_walker( imInput, seedIndices, seedLabels );

timeElapsed = toc

segMaskRGB = label2rgbND(segMask);

% display result
objSeedMaskRGB = label2rgbND(objSeedMask);
imseriesmaskshowrgb(imInput, {segMaskRGB, objSeedMaskRGB} );
set( gcf, 'Name', 'Random Walker Segmentaion Result -- mask, object-seed-mask' );

if ndims(imInput) == 2
    
    segMaskRep = label2rgbND( repmat(segMask, [ones(1,ndims(segMask)), numLabels]) );
    objSeedMaskRep = label2rgbND( repmat(objSeedMask, [ones(1,ndims(objSeedMask)), numLabels]) );

    imseriesmaskshowrgb(pixToLabelProbabilityMap * 100, {segMaskRep, objSeedMaskRep} ); 
    set( gcf, 'Name', 'Random Walker Segmentaion Result -- pixelLabelProbabilityMap' );
    
end

