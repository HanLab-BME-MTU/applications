function [segMask, pixToLabelProbabilityMap] = random_walker_segmentation( imInput, seedIndices, seedLabels, varargin )
% random_walker_segmentation: An ND implementation of the random walker
% segmentation algorithm proposed by Grady et. al. (PAMI 2006].
%
%     [segMask, pixToLabelProbabilityMap, labelIDList] = random_walker_segmentation( imInput, seedIndices, seedLabels, varargin ); 
%
%     Required Input Arguments:
% 
%                         imInput: input ND image
% 
%                     seedIndices: linear pixel indices (from sub2ind) of the 
%                                  seed pixels
% 
%                      seedLabels: label ids of the seed pixels specified in 
%                                  seedIndices argument. Both seedIndices and 
%                                  seedLabels must have the same size.
% 
%     Optional Input Arguments:
% 
% 
%                         spacing: pixel spacing of the input image.
%                                  can be a scalar value or an array of size
%                                  equal to the number of image dimensions.
%                                  Default: 1 (isotropic spacing)                      
% 
%                            beta: scalar value used in the gaussian weighting 
%                                  function used to compute the edge weights.
%                                  Default: 90
% 
%                       neighconn: specifies the connectivity type of the 
%                                  neighborhood graph.
%                                  Default: 0
%                                  Options: 0 or 1
% 
%                                   0 - 4-connected for 2D image, 
%                                       6-connected for 3D image,
%                                       2N-connected for ND image
% 
%                                   1 - 8-connected for 2D image,
%                                       26-connected for 3D image,
%                                       (3^N-1)-connected for ND image                              
% 
%     Output Arguments:
% 
%                        segMask: output segmentation mask. 
%                                 pixels are assigned to the label with highest 
%                                 probability.
% 
%       pixToLabelProbabilityMap: probabilities of pixels belonging to each
%                                 segmentation label.
%                      
% 
% Example:  
% 
%        % generate image
%        imInput = zeros( 100 , 100 );
%        
%        gmobj_bgnd_ex = gmdistribution( 50 , 30 );
%        imInput_bgnd_ex = reshape( random( gmobj_bgnd_ex , numel( imInput ) ) , size( imInput ) );
% 
%        gmobj_fgnd_ex = gmdistribution( 90 , 30 );
%        imInput_fgnd_ex = reshape( random( gmobj_fgnd_ex , numel( imInput ) ) , size( imInput ) );
%        
%        imInput = imInput_bgnd_ex;
%        
%        image_size = size( imInput );
%        
%        imInput( image_size(1) * 0.25 : image_size(1) * 0.75 , image_size(2) * 0.25 : image_size(2) * 0.75 ) = imInput_fgnd_ex( image_size(1) * 0.25 : image_size(1) * 0.75 , image_size(2) * 0.25 : image_size(2) * 0.75 );
%        
%        imInput = round( imInput );
% 
%        % get seed points
%        [ fgnd_seed_points , bgnd_seed_points ] = get_fgnd_bgnd_seeds_2d_strokes( imInput, DisplayRange );
%        
%        fgndSeedIndices = sub2ind( size(imInput), fliplr(fgnd_seed_points) );
%        fgndSeedMask = zeros(size(imInput));
%        fgndSeedMask(fgndSeedIndices) = 1;
% 
%        bgndSeedIndices = sub2ind( size(imInput), fliplr(bgnd_seed_points) );
%        bgndSeedMask = zeros(size(imInput));
%        bgndSeedMask(bgndSeedIndices) = 1;
%        
%        seedIndices = [ fgndSeedIndices; bgndSeedIndices ];
%        seedLabels = [ones(size(fgndSeedIndices)); zeros(size(bgndSeedIndices))];
%        
%        % run random walker segmentation algorithm
%        [segMask, pixToLabelProbabilityMap] = random_walker_segmentation( imInput, seedIndices, seedLabels );
%        
%        % display result
%        imseriesmaskshow(imInput, {segMask, fgndSeedMask, bgndSeedMask}, 'maskColors', [0 1 0; 0 0 1; 1 0 0] );
%        set( gcf, 'Name', 'Random Walker Segmentaion Result -- mask, fgndSeeds, bgndSeeds' );
%        
%        imseriesmaskshow(pixToLabelProbabilityMap, {segMask, fgndSeedMask, bgndSeedMask}, 'maskColors', [0 1 0; 0 0 1; 1 0 0] );
%        set( gcf, 'Name', 'Random Walker Segmentaion Result -- pixelLabelProbabilityMap, fgndSeeds, bgndSeeds' );
% 
% Author: Deepak Roy Chittajallu
% 
% References:
% 
% [1] Grady, L. (2006). "Random Walks for Image Segmentation." IEEE Transactions on Pattern Analysis and Machine Intelligence, 28(11): 1768-1783.
% 

    p = inputParser;
    p.CaseSensitive = false;
    p.addRequired( 'imInput', @(x) (ndims(x) > 1) );
    p.parse( imInput );
    
    imdims = ndims(imInput);   
    
    p.addParamValue( 'spacing', ones(1, imdims), @(x) (isnumeric(x) && ismember(numel(x), [1,imdims])) );
    p.addParamValue( 'sigma', 0.5 / sqrt(90), @isscalar );
    p.addParamValue( 'neighconn', 0, @(x) (isscalar(x) && ismember(x, [0,1])));
    p.addParamValue( 'flagNormalize', true, @(x) (isscalar(x) && islogical(x)));
    p.addParamValue( 'flagParallelize', true, @(x) (isscalar(x) && islogical(x)));
    p.parse(imInput, varargin{:} );
    
    parameters = p.Results;
    
    % image size
    imsize = size(imInput);
        
    % build grid graph with the specified neighborhood connectivity
    [~, edges, edgeType, edgeTypeOffset] = buildGridGraphND( imsize, parameters.neighconn );

    % assign weights to edges
    I1 = imInput(edges(:,1));
    I2 = imInput(edges(:,2));
    %edgeIntensityDiff = abs(I1 - I2);
    edgeIntensityDiff = abs(I1 - I2) ./ (eps + min([I1, I2], [], 2));

    if parameters.flagNormalize
        edgeIntensityDiff = mat2gray(edgeIntensityDiff);
    end

    spacingNmzed = parameters.spacing / max(parameters.spacing);
    edgeOffsetDist = edgeTypeOffset * diag(spacingNmzed);
    edgeOffsetDist = sqrt(sum(edgeOffsetDist.^2, 2));
    edgeDist = edgeOffsetDist(edgeType);
        
    edgeWeights = exp( -(2.0*parameters.sigma^2)^-1 * edgeIntensityDiff.^2 ./ edgeDist.^2 ) + 1e-5;

    % construct laplacian matrix
    adjMatrix = sparse( [edges(:,1); edges(:,2)], ...
                        [edges(:,2); edges(:,1)], ...
                        [edgeWeights; edgeWeights], ...
                        numel(imInput), numel(imInput) );
                    
    L = diag(sum(adjMatrix)) - adjMatrix;
    
    % setup the dirichlet problem
    labelIDList = unique(seedLabels);
    numLabels = numel(labelIDList);    
    boundaryCondition = zeros(numel(seedIndices), numLabels);
    
    for k = 1:numLabels
       boundaryCondition(seedLabels == labelIDList(k), k) = 1;
    end
   
    % solve the dirichlet problem    
    unlabeledPixelIndices = (1:numel(imInput))';
    unlabeledPixelIndices(seedIndices) = [];
    B = -L(unlabeledPixelIndices, seedIndices) * boundaryCondition;
    A = L(unlabeledPixelIndices, unlabeledPixelIndices);
    
    %spparms ('spumoni',3)
    
    switch imdims
    
        case 2
            
            X = A \ B;
            
        case 3
            
            M = diag(diag(A));
            X = zeros(numel(unlabeledPixelIndices), numLabels);
            
            maxiter = 5000;
            tol = 1e-5;
            
            if parameters.flagParallelize
                
                flagPoolOpenedAlready = matlabpool( 'size' ) > 0;        
                if ~flagPoolOpenedAlready 
                    matlabpool open;
                end            

                parfor i = 1:numLabels                                        
                    X(:,i) = pcg(A, B(:,i), tol, maxiter, M);  
                end
                
                if parameters.flagParallelize && ~flagPoolOpenedAlready
                    matlabpool close;
                end    
                
            else
                
                for i = 1:numLabels                                        
                    X(:,i) = pcg(A, B(:,i), tol, maxiter, M);  
                end
                
            end

    end
    
    %spparms ('spumoni',0)
    
    pixToLabelProbabilityMap = zeros(numel(imInput), numLabels);
    pixToLabelProbabilityMap(seedIndices, :) = boundaryCondition;
    pixToLabelProbabilityMap(unlabeledPixelIndices, :) = X;
    
    % generate segmentation mask by assigning pixels to the label with highest probability
    [maxprob, segMask] = max(pixToLabelProbabilityMap, [], 2);
    segMask = labelIDList(segMask);
    
    % reshape stuff
    segMask = reshape(segMask, imsize);
    pixToLabelProbabilityMap = reshape(pixToLabelProbabilityMap, [imsize, numLabels]);
    
end


