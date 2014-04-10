function [varargout] = GVFND(imFeature, weightReg, numIterations, varargin )
% Computes the Gradient Vector Field in ND
%

    p = inputParser;
    p.addRequired( 'imFeature', @(x) (isnumeric(x) && ndims(x) >= 2) );
    p.addRequired( 'weigthReg', @(x) (isnumeric(x) && isscalar(x)) );
    p.addRequired( 'numIterations', @(x) (isnumeric(x) && isscalar(x)) );
    p.parse(imFeature, weightReg, numIterations );
    
    p.addParamValue('spacing', ones(1, ndims(imFeature)), @(x) ( isnumeric(x) && ~isscalar(x) && numel(x) == ndims(imFeature) ) );
    p.addParamValue('flagDebugMode', false, @(x) ( islogical(x) & isscalar(x) ) );
    p.addParamValue('flagParallelMode', false, @(x) ( islogical(x) & isscalar(x) ) );
    p.parse(imFeature, weightReg, numIterations, varargin{:} );
    
    spacing = p.Results.spacing;
    flagDebugMode = p.Results.flagDebugMode;
    flagParallelize = p.Results.flagParallelMode;
    
    imdims = ndims(imFeature);
    imsize = size(imFeature);       
    spacing = num2cell(spacing);
    
    if flagParallelize
        flagPoolOpenedAlready = matlabpool( 'size' ) > 0;        
        if ~flagPoolOpenedAlready 
            matlabpool open;
        end            
    end    
    
    % pad the feature image so that pixels on the boundary are handled
    % properly
    padMethod = 'symmetric';
    padSize = ones(1,imdims); % just enough to compute the laplacian at the boundary
    imFeaturePadded = padarray(imFeature, padSize,padMethod);
    unpadmask = padarray(ones(imsize), padSize, 0 );
    borderCorrectionInd = getMatrixPaddingIndices(imsize, padSize, padMethod, 'both' );
        
    for j = 1:imdims
        borderCorrectionInd{j} = borderCorrectionInd{j} + 1;
    end
    
    % normalize the feature image to the range [0,1]
    imFeaturePadded = (imFeaturePadded - min(imFeaturePadded(:))) / (max(imFeaturePadded(:)) - min(imFeaturePadded(:)));
    
    % initialize vector field
    vectorField = cell(1, imdims);    
    [vectorField{:}] = gradient(imFeaturePadded, spacing{:});
    initGVF = vectorField;
    
    gmag2 = zeros(size(imFeaturePadded));
    for j = 1:imdims
        gmag2 = gmag2 + vectorField{j} .* vectorField{j};
    end
    
    % Iteratively diffuse the vector field
    fprintf(1, '\n' );
    iterChange = zeros(1,imdims);   
    for i = 1:numIterations        
        if flagParallelize
            parfor j = 1:imdims

                % boundary correction/ensure-mirroring
                vold = vectorField{j}(borderCorrectionInd{:}); 

                % compute euler step
                vectorField{j} = vold + weightReg * (2 * imdims) * del2(vectorField{j}, spacing{:}) - gmag2 .* (vectorField{j} - initGVF{j});

                % record maximum change from previous iteration
                vchange = vectorField{j} - vold;            
                iterChange(j) = max(abs(vchange(:)));            

            end                                
        else            
            for j = 1:imdims

                % boundary correction/ensure-mirroring
                vold = vectorField{j}(borderCorrectionInd{:}); 

                % compute euler step
                vectorField{j} = vold + weightReg * (2 * imdims) * del2(vectorField{j}, spacing{:}) - gmag2 .* (vectorField{j} - initGVF{j});

                % record maximum change from previous iteration
                vchange = vectorField{j} - vold;            
                iterChange(j) = max(abs(vchange(:)));            

            end        
        end
        % print iteration progress
        fprintf(1, '%3d ', i );
        if( rem(i,20) == 0 )
            fprintf(1, '\n' );
        end                        
        
        changeTrend(i) = min(iterChange);
        
        if ~any(iterChange > 10^-4)
            fprintf(1, '\nsteady-state reached in %d iterations\n', i );
            break;
        end
    end    
    fprintf(1, '\n' );
    
    if flagDebugMode
        figure, plot(changeTrend);        
    end
    
    % unpad the components of the vector field
    for j = 1:imdims
        vectorField{j} = reshape( vectorField{j}(unpadmask > 0), imsize );
    end    
    
    varargout = vectorField;
    
    if flagParallelize && ~flagPoolOpenedAlready
        matlabpool close;
    end
    
end