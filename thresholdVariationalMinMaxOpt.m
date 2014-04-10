function [ imThresholdSurface, varargout ] = thresholdVariationalMinMaxOpt( imInput, varargin )
% Computes a threshold using variational minimix optimization.
% 
% [ imThresholdSurface ] = thresholdVariationalMinMaxOpt( imInput, varargin )
% [ imThresholdSurface, imMask ] = thresholdVariationalMinMaxOpt( imInput, varargin )
% 
% Please refer to the following paper on further details about how this
% algorithm works:
% 
% Saha, B. N. and N. Ray (2009). "Image thresholding by variational minimax 
% optimization." Pattern Recognition 42(5): 843-856.
% 
%     Required Input Arguments:     
%         
%                          imInput: Input ND Image
%                          
%     Optional Input Arguments:    
% 
%                          spacing: pixel spacing of the input image.
%                                   
%                                   must be an array of size equal to the 
%                                   number of image dimensions.
%                                   Default: isotropic spacing is assumed
%                                   
%     Optional Param/Value Arguments:    
%                                
%                                q: an exponential weight of the gradient 
%                                   term h in the paper.
%                                   
%                                   Must be a positive scalar value.
%                                   Default: 1.0
%                                   
%                        tolerance: specifies a convergence tolerance for 
%                                   the energy being minimized.
%                                   
%                                   While minimizing the energy using gradient 
%                                   descent procedure, if the change in energy 
%                                   from the last iteration is below this value 
%                                   then convergence assumed and gradient
%                                   descent is stopped.
% 
%                                   Must be a positive scalar value.
%                                   Default: 1e-5
%                                   
%                    maxIterations: specifies the maximum number of iterations             
%                                   for the gradient descent procedure.
%                                       
%                                   To find a good value for this, run this
%                                   function in debug mode on one of the 
%                                   images from your database and look at 
%                                   the displayed energy vs iteration plot.
%                                   And pick a value when it seems to 
%                                   converge on your images.
%                                   
%                                   Must be an positve integer.
%                                   Default: 1000
%                                   
%                                   
%                        debugMode: true/false
%                                   A bunch of stuff is printed and displayed
%                                   in debug mode  
%                                   Default: false
%                            
%     Output Arguments:
% 
%               imThresholdSurface: A matrix of the same size as the input 
%                                   image where in each element contains the 
%                                   threshold value at that pixel.
% 
%                imMask (optional): A binary foreground mask obtained by 
%                                   applying the computed threshold to the 
%                                   input image
%                                   
% Examples:
% 
%     imInput = imread('coins.png');
%     
%     [imThresholdSurface, imMask] = thresholdVariationalMinMaxOpt( imInput, 'debugMode', true );
%   
%     
% Author: Deepak roy Chittajallu (Created Mar 11, 2013)
% 
                                  
    p = inputParser;
    p.addRequired( 'imInput', @(x) (ismember(ndims(x), [2,3])) );
    p.parse( imInput );
    
    numDims = ndims(imInput);
    
    p.addOptional('spacing', ones(1,numDims), @(x) (numel(x) == numDims));    
    p.addParamValue('q', 1.0, @isscalar); % same as q in paper
    p.addParamValue('tolerance', 1e-10, @isscalar );
    p.addParamValue('maxIterations', 1000, @isscalar);        
    p.addParamValue('debugMode', true, @islogical);    
    p.parse(imInput, varargin{:} );

    funcParameters = p.Results;

    if funcParameters.debugMode
        funcParameters
    end
    
    imInput = double(imInput);
    inputIntensityRange = [min(imInput(:)), max(imInput(:))];
    cellSpacing = num2cell(funcParameters.spacing);
    
    % Standardize intensity range
    imInput = mat2gray(imInput);
    
    % Compute gradient
    delI = delMat(imInput, funcParameters.spacing);
    magDelI = sqrt(delMagSq(delI));
    
    % Pre-compute H - strength of gradient at each pixel/voxel
    imH = magDelI.^funcParameters.q;
    imH = imH / max(imH(:));
    
    % Compute the threshold surface
    otsuThresh = thresholdOtsu(imInput);    
    imT = zeros(size(imInput)) + otsuThresh;       
    
    if funcParameters.debugMode
        fprintf( '\nComputing threshold using variational minimax optimization ...\n');
        hVis = figure;
        if numDims == 2
            imOtsuMask = imInput > otsuThresh;
            subplot(2,3,1), imshow(imInput, []); title( 'Input Image' );                    
            subplot(2,3,2), imshow(imOtsuMask, []); title( 'Otsu Result' );
            subplot(2,3,3)
                title('Energy vs Iterations');
                xlabel('Iterations');
                ylabel('Energy');
        end
    end
    
    for i = 1:funcParameters.maxIterations

        % Compute term E1 in the cost function
        E1 = Compute_E1(imInput, imT, imH);     
        
        % Compute term E2 in the cost function - regularizes the threshold 
        [E2, delT] = Compute_E2(imT, funcParameters.spacing);
        
        % Compute optimal alpha
        alpha = E2 / norm([E1, E2]);
        
        % Compute Total Energy
        w1 = sqrt(1 - alpha^2);
        w2 = alpha;
        E = w1 * E1 + w2 * E2;

        if funcParameters.debugMode
            fprintf( '\nIteration %d/%d: E = %f, alpha = %f', ...
                     i, funcParameters.maxIterations, E, alpha);           
                 
            if numDims == 2
                
                figure(hVis);
                
                subplot(2,3,3)
                    hold on;
                        if i > 1
                           plot([i-1; i], [oldE; E], 'b-', 'LineWidth', 2.0); 
                           drawnow; 
                        end
                    hold off;                
                
                subplot(2,3,4)
                    imshow(imT, []);
                    title(sprintf('Threshold Surface - Iteration %d', i));
                    xlabel(sprintf( 'E = %.3f, alpha = %.3f', E, alpha));
                    
                subplot(2,3,5), imshow(imInput > imT, []);    
                    title('Foreground Mask');
                    
                subplot(2,3,6)
                    imDiffWithOtsu = (imInput > imT) - imOtsuMask;
                    imshow(imDiffWithOtsu, []);    
                    title('Difference with Otsu');
                    
                drawnow;
                
            else
                figure(hVis);
                hold on;
                    if i > 2
                       plot([i-1; i], [oldE; E], 'b-', 'LineWidth', 2.0); 
                       drawnow; 
                    end
                hold off;                
            end
        end
        
        if i > 1 && abs(alpha - oldAlpha) < funcParameters.tolerance
            fprintf( '\nConvergence reached after %d iterations !!!\n', i);
            break;
        else
            oldAlpha = alpha;
            oldE = E;
        end
        
        % update threshold by stepping along the steepest descent direction
        
            % Compute direction of steepest descent - obtained by solving
            % the euler-lagrange equation
            lapT = (2 * numDims) * del2(imT, cellSpacing{:});
            dirT = (w1 * imH .* (imInput - imT)) + (w2 * lapT);
            
            % Compute the optimal step size using line search
            delDirT = delMat(dirT, funcParameters.spacing);
            
            numVal = w1 * sumall(imH .* (imInput - imT) .* dirT) - w2 * sumall(dotDel(delT, delDirT));
            denomVal = w1 * sumall(imH .* dirT.^2) + w2 * sumall(delMagSq(delDirT));            
            stepSize = numVal / denomVal;
            
            % update threshold
            imT = imT + (stepSize * dirT);     
            
        if funcParameters.debugMode
            fprintf(', step-size = %f', stepSize);
        end
        
    end
    
    if i == funcParameters.maxIterations
        fprintf( '\nMaximum iteration value of %d has been reached before convergence\n', ...
                 funcParameters.maxIterations);
    end
    
    imThresholdSurface = inputIntensityRange(1) ...
                         + imT * (inputIntensityRange(2) - inputIntensityRange(1));    
    
    if nargout > 1
        imMask = double(imInput > imT);
        varargout{1} = imMask;        
    end
end

% Compute term E1 in cost function - quantifies the error between the image
% and the threshold surface
function [E1] = Compute_E1( I, T, H )

    E1 = 0.5 * sumall(H .* (I - T).^2);
    
end

% Compute term E2 in cost function - regularizes the threshold surface
% and the threshold surface
function [E2, varargout] = Compute_E2( T, spacing )

    numDims = ndims(T);    
    if ~exist( 'spacing', 'var' )
        spacing = ones(1,numDims);
    end
    
    % Compute del(T)
    delT = delMat(T, spacing);
    
    % Compute magnitude of del(T)
    magSqDelT = delMagSq(delT);
    
    % Compute E2
    E2 = 0.5 * sumall(magSqDelT);

    % return optional output arguments if requested
    if nargout > 1
        varargout{1} = delT;
    end
    
    if nargout > 2
        varargout{2} = magSqDelT;
    end
    
end

% Sum all elements in a matrix
function [M_sum] = sumall( M )
    M_sum = sum(M(:));
end

% Find dot product between grad of two ND matrices
function [dotDel_val] = dotDel(delMat1, delMat2)

    numDims = ndims(delMat1);
    
    dotDel_val = 0;
    for i = 1:numDims
        dotDel_val = dotDel_val + delMat1{i} .* delMat2{i};
    end
    
end

% Compute del of a matrix
function [del_val] = delMat(M, spacing)

    numDims = ndims(M);
    if ~exist( 'spacing', 'var' )
        spacing = ones(1,numDims);
    end

    cellSpacing = num2cell( spacing );
    
    % compute gradient along each direction
    del_val = cell(1,numDims);
    [del_val{:}] = gradient(M, cellSpacing{:});    
    
%     % zero-out borders
%     cellAllElemInd = cell(1,numDims);
%     for i = 1:numDims
%         cellAllElemInd{i} = ':';
%     end
%     
%     graddir = [2, 1, 3:numDims];
%     for i = 1:numDims
%         
%         % clear first
%         cellBorderInd = cellAllElemInd;
%         cellBorderInd{graddir(i)} = 1;
%         del_val{i}(cellBorderInd{:}) = 0;
%         
%         % clear last
%         cellBorderInd = cellAllElemInd;
%         cellBorderInd{graddir(i)} = size(M, graddir(i));
%         del_val{i}(cellBorderInd{:}) = 0;
%         
%     end
    
end

% Compute the squared magnitude of the del of a matrix
function [magSq] = delMagSq( delM )

    magSq = dotDel(delM, delM);
    
end
