function [blobStats, imBlobSeedMask, imBlobSegMask, PARAMETERS] = segmentFociInsideNuclei(imInput, blobDiameterRange, varargin)
%% Detect blobs of varying size using a SIFT like approach
%
% Author: Deepak Roy Chittajallu
%

%% PARSE PARAMETERS

    p = inputParser;
    p.addRequired( 'imInput', @(x) (isnumeric(x) && ismember( ndims(x), [2,3] )) );
    p.addRequired( 'blobDiameterRange', @(x) (numel(x) == 2) );  
    p.addParamValue( 'spacing', ones( 1, ndims(imInput) ), @(x) (isnumeric(x) && numel(x) == ndims(imInput)) );
    
    p.addParamValue( 'sigmaLogDelta', 0.2447, @(x) @(x)(isnumeric(x) && isscalar(x)) );
    p.addParamValue( 'minSignalToBackgroundRatio', 1.0, @(x) (isnumeric(x) && isscalar(x) && x >= 1) );
    p.addParamValue( 'minDoGResponse', eps, @(x) (isnumeric(x) && isscalar(x) && x >= 0) );
    p.addParamValue( 'minHessianEigenRatio21', 0.0, @(x) (isnumeric(x) && isscalar(x) && x >= 0 && x <= 1.0) );
    p.addParamValue( 'minHessianEigenRatio31', 0.0, @(x) (isnumeric(x) && isscalar(x) && x >= 0 && x <= 1.0) );
    
    p.addParamValue( 'roiMask', [], @(x) ( (isnumeric(x) || islogical(x)) && ndims(x) == ndims(imInput) && all(size(x) == size(imInput)) ) );    
    p.addParamValue( 'minDistanceToROIBoundary', 1.0, @(x)(isnumeric(x) && isscalar(x) && x >= 1.0) );
    p.addParamValue( 'fociDetectionModelFile', [], @(x) (isempty(x) || (ischar(x) && exist(x, 'file'))) );
    p.parse(imInput, blobDiameterRange, varargin{:});
    
    PARAMETERS = p.Results;

%% Good parameters emipirically found for different tasks
% Nuclei seed point detection
%     blobDiameterRange = [8, 20];
%     PARAMETERS.numScales = 10;
%     PARAMETERS.minDoGResponse = 0.003;
%     PARAMETERS.minSignalToBackgroundRatio = 3.0;
%     PARAMETERS.minHessianEigenRatio21 = 0.5;
%     PARAMETERS.minHessianEigenRatio31 = 0.0;
%     PARAMETERS.minDistanceToROIBoundary = 0.0;  
%      
% Foci detection inside nuclei
%     blobDiameterRange = min(PARAMETERS.spacing) * [3, 7];
%     PARAMETERS.numScales = 5;
%     PARAMETERS.minDoGResponse = 1e-3;
%     PARAMETERS.minSignalToBackgroundRatio = 1.5;
%     % maxHessianEigenRatio12 - med = 1.44, mad = 2, sigma = 0.7420 -> s1 = 2.15, s2 = 2.85, s3 = 3.55, s4 = 4.41, s5 = 5.15
%     PARAMETERS.minHessianEigenRatio21 = 0.5; 
%     PARAMETERS.minHessianEigenRatio31 = 0.00;
%     PARAMETERS.minDistanceToROIBoundary = 1.0;
%
%%

    fociDetectionTimer = tic;
    imdims = ndims(imInput);    
    
    % Pre-processing
    fprintf('\n>> Preprocessing ... \n');

    imPreprocessed = mat2gray(imInput); % standardize
%     imPreprocessed = matitk('FMEDIAN', [1, 1, zeros(1,imdims-2)], imPreprocessed); % applying in each plane because of intensity attenuation with depth
    imPreprocessed = filterGaussND(imPreprocessed, 0.5 * min(PARAMETERS.spacing), ...
                                   'spacing', PARAMETERS.spacing);
    
    if ~isempty(PARAMETERS.roiMask)
        distToROIBoundary = bwdistsc(~PARAMETERS.roiMask, PARAMETERS.spacing);
    end
        
    % Enhance foci usint a tophat filter
    fprintf('\n>> Estimating local background to compute the contrast of blobs against background ... \n');

    krnlMax = streldisknd( round(0.5 * max(blobDiameterRange) ./ PARAMETERS.spacing(1:2)) ); 
    imLocalBackground = imopen(imPreprocessed, krnlMax);
    imSignalToBackgroundRatio = imPreprocessed ./ (eps + imLocalBackground);

    % detect blobs in scale space
    fprintf('\n>> Detecting blobs at multiple scales ... \n');

    blobLocIndices = [];
    blobRadii = [];
    blobEigenH = [];
    blobLoGResponse = [];

    sigmaLogDelta = PARAMETERS.sigmaLogDelta;
    sigmaLogRange = log2( (0.5 * sort(blobDiameterRange) / sqrt(imdims)) );
    numScales = ceil( (sigmaLogRange(2) - sigmaLogRange(1)) / sigmaLogDelta );
    sigmaInit = 2^(sigmaLogRange(1)-sigmaLogDelta);
    sigmaFunc = @(s) (2.^(s*sigmaLogDelta) * sigmaInit );
    c = sqrt( 2^(2*sigmaLogDelta) - 1 );

    imG_prev = filterGaussND(imPreprocessed, sigmaInit, 'spacing', PARAMETERS.spacing); 
    dogScaleWindow = cell(1,2);

    for s = 1:(numScales+2)

        curSigmaStep = c * sigmaFunc(s-1); 
        imG = filterGaussND(imG_prev, curSigmaStep, 'spacing', PARAMETERS.spacing);

        % Compute LoG using DoG approximation
        imDoG_prev = imG - imG_prev; 

        % Compute local maxima for scale (s-2)
        if s <= 2
            dogScaleWindow{s} = imDoG_prev;
        else                

            dogScaleWindow{3} = imDoG_prev;    
            rBlob = sigmaFunc(s-2) * sqrt(imdims);
            fprintf('\nDetecting blobs of diameter %.2f um ... ', 2.0 * rBlob);
            curFociDetectionTimer = tic;

            % compute local maxima in scale space
            imMultiscaleLoG = max(-cat(imdims+1, dogScaleWindow{:}), [], imdims+1);
            w = 2 * ceil(rBlob ./ PARAMETERS.spacing) + 1;
            if imdims == 3
                imLocalMax = locmax3d(imMultiscaleLoG, w, 'ClearBorder', false);           
            else
                imLocalMax = locmax2d(imMultiscaleLoG, w, 1);           
            end
            imLocalMax(imLocalMax ~= -dogScaleWindow{2}) = 0;

            curBlobLocIndices = find(imLocalMax > 0);
            
            % prune blobs outside a specified roi mask
            if ~isempty(PARAMETERS.roiMask)
                curBlobLocIndices(~PARAMETERS.roiMask(curBlobLocIndices)) = [];
            end
            
            % compute eigen values of the hessian matrix
            if ~isempty(curBlobLocIndices)

                % compute hessian
                imCur = padarrayXT(imG_prev - dogScaleWindow{2}, 2 * ones(1, imdims), 'symmetric');
                x = 2 + ind2submat(size(imLocalMax), curBlobLocIndices);
                H = hessianAt(imCur, x, PARAMETERS.spacing);

                % compute eigen values of the hessian
                numPoints = numel(curBlobLocIndices);                        
                curBlobEigenH = zeros(numPoints, imdims);
                for i = 1:numPoints
                    curBlobEigenH(i, :) = eig( H(:,:,i) );
                end

            end            

            % update blob list
            if numel(curBlobLocIndices) > 0
                curBlobRadii = repmat(rBlob, size(curBlobLocIndices));
                blobLocIndices = cat(1, blobLocIndices, curBlobLocIndices);
                blobRadii = cat(1, blobRadii, curBlobRadii);
                blobEigenH = cat(1, blobEigenH, curBlobEigenH);
                blobLoGResponse = cat(1, blobLoGResponse, imLocalMax(curBlobLocIndices));
            end
            
            % update DoG scale window
            dogScaleWindow(1) = [];

            % time elapsed
            fprintf( 'took %.2f seconds\n', toc(curFociDetectionTimer));
        end

        imG_prev = imG;

    end

    % Blob stats
    blobStats = [];

    if ~isempty(blobLocIndices)

        blobStats = struct('PixelLocationIndex', num2cell(blobLocIndices), ...
                           'Radius', num2cell(blobRadii), ...
                           'SignalToBackgroundRatio', num2cell(imSignalToBackgroundRatio(blobLocIndices)), ...
                           'DoGResponse', num2cell(blobLoGResponse));

        for i = 1:numel(blobStats)
            blobStats(i).HessianEigenValues = blobEigenH(i,:);
        end
        
        if ~isempty(PARAMETERS.roiMask)
            blobDist = num2cell(distToROIBoundary(blobLocIndices));
            [blobStats.DistanceToBoundary] = blobDist{:};
        end

    end

    % prune blobs based on various criterion
    if ~isempty(blobLocIndices)

        flagValidBlobs = true(size(blobLocIndices));
        
        % Apply cutoffs
        fprintf('\n>> Pruning invalid foci using the specified cutoffs ... \n');

        % prune blobs where input signal is below a specified cutoff
        flagValidBlobs(imSignalToBackgroundRatio(blobLocIndices) < PARAMETERS.minSignalToBackgroundRatio) = false;

        % prune blobs with LoG response below a specified cutoff
        flagValidBlobs(blobLoGResponse < PARAMETERS.minDoGResponse) = false;

        % prune edge responses by analyzing the eigen-values of the hessian matrix
        blobEigenRatio21 = blobEigenH(:,2) ./ (eps + blobEigenH(:,1));
        blobEigenRatio31 = blobEigenH(:,3) ./ (eps + blobEigenH(:,1));
        flagValidBlobs(sum(sign(blobEigenH), 2) ~= -3 | ...
                       blobEigenRatio21 <= PARAMETERS.minHessianEigenRatio21 | ...
                       blobEigenRatio31 <= PARAMETERS.minHessianEigenRatio31) = false; 
        
        % prune blobs too close to roi boundary
        % This gets rid of quite a few edge responses but it
        % assumes that the roi boundary coincides with edges in the image
        if ~isempty(PARAMETERS.roiMask) && PARAMETERS.minDistanceToROIBoundary > 0
            flagValidBlobs(distToROIBoundary(blobLocIndices) < PARAMETERS.minDistanceToROIBoundary * blobRadii) = false;
        end
                   
        % prune
        blobStats(~flagValidBlobs) = [];
        blobLocIndices(~flagValidBlobs) = [];
        
    end
    
    % apply foci detection model if provided
    if ~isempty(PARAMETERS.fociDetectionModelFile)

       fprintf('\n>> Pruning invalid foci using the specified learning model ... \n');

       flagValidBlobs = true(size(blobLocIndices));
       
       % weka imports
       import weka.core.*;    
       import weka.core.SerializationHelper.*;

       wekaModel = weka.core.SerializationHelper.readAll( PARAMETERS.fociDetectionModelFile );

       fprintf( '\nProgress: \n' );

       last_percent_done = 0;
       numPrint = 0;

        for bid = 1:numel(blobStats)

            featureStruct = ComputeFociDetectionFeatures(imInput, blobStats, bid, ...
                                                         'spacing', PARAMETERS.spacing);                                                         
            [featureVec , featureNameList] = ConvertFeatureStructToFeatureVec( featureStruct );                                                 

            predictedLabel = wekaPredictInstance(wekaModel, featureVec, featureNameList);

            if strcmpi(predictedLabel, 'Bad_Detection')
                flagValidBlobs(bid) = false;
            end 

            percent_done = round(100 * bid / numel(blobStats));       

            if percent_done > last_percent_done
                fprintf( '%.2d%%  ', percent_done );
                last_percent_done = percent_done;
                numPrint = numPrint + 1;
                if mod( numPrint, 10 ) == 0
                   fprintf( '\n' ); 
                end
            end        

        end
            
        % prune
        blobStats(~flagValidBlobs) = [];
        blobLocIndices(~flagValidBlobs) = [];
        
    end
    
    % Build blob seed mask
    if nargout > 1
        imBlobSeedMask  = zeros( size(imInput) );
        imBlobSeedMask(blobLocIndices) = 1:numel(blobLocIndices);
    end

    % Build blob segmentation mask based on size
    if nargout > 2
        seedPos = ind2submat(size(imInput), blobLocIndices) * diag(PARAMETERS.spacing);
        kd = KDTreeSearcher(seedPos);

        pixelPos = ind2submat(size(imInput), (1:numel(imInput))') * diag(PARAMETERS.spacing);
        [closestSeedInd, distanceToSeed] = kd.knnsearch(pixelPos);

        imBlobSegMask = zeros( size(imInput) );
        flagIsPixelInSeedVicinity = distanceToSeed <= blobRadii(closestSeedInd);
        imBlobSegMask( flagIsPixelInSeedVicinity ) = closestSeedInd( flagIsPixelInSeedVicinity );
        imBlobSegMask = imdilate(imBlobSegMask, streldisknd(ones(1,imdims)));                
    end

    % Done
    fprintf('\nblob detection took a total of %.2f seconds\n', toc(fociDetectionTimer) );
    
    fprintf('\nA total of %d blobs were detected\n', numel(blobStats));
        
end

function [H] = hessianAt(f, points, spacing)
%
%   REF: https://v8doc.sas.com/sashtml/ormp/chap5/sect28.htm
%
    if ~exist('spacing', 'var')
        spacing = ones(1, ndims(f));
    end

    fdims = ndims(f);    
    numPoints = size(points,1);
    
    H = zeros([fdims, fdims, numPoints]);

    pointsPlus = cell(1, fdims);
    pointsMinus = cell(1, fdims);
    for i = 1:fdims
        pointsPlus{i} = points;
        pointsPlus{i}(:,i) = points(:,i) + 1;
        
        pointsMinus{i} = points;
        pointsMinus{i}(:,i) = points(:,i) - 1;        
    end
    
    g = @(p, i) gradAt(f, i, p, spacing);
    
    for i = 1:fdims
        for j = i:fdims
            H(i,j,:) = ((g(pointsPlus{j}, i) - g(pointsMinus{j}, i)) / (4 * spacing(j))) + ...
                       ((g(pointsPlus{i}, j) - g(pointsMinus{i}, j)) / (4 * spacing(i)));
            if i ~= j
                H(j,i,:) = H(i,j,:);
            end            
        end
    end
    
end

function gf = gradAt(f, dir, points, spacing)
%
%   REF: https://v8doc.sas.com/sashtml/ormp/chap5/sect28.htm
%

    if ~exist('spacing', 'var')
        spacing = ones(1, ndims(f));
    end

    pointsPlus = points;
    pointsPlus(:, dir) = points(:, dir) + 1;
    pointsPlus = submat2ind(size(f), pointsPlus);
    
    pointsMinus = points;
    pointsMinus(:, dir) = points(:, dir) - 1;
    pointsMinus = submat2ind(size(f), pointsMinus);
    
    gf = (f(pointsPlus) - f(pointsMinus)) / (2 * spacing(dir));
    
end
