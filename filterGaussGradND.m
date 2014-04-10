function [ imGaussGrad ] = filterGaussGradND(im, sigma, direction, varargin )
% filterGaussGradND: An ND implementation of the Derivative of Gaussian Filter
% 
%     [ imGaussGrad ] = filterGaussGradND(im, sigma, direction, varargin );
% 
%     Required Input Arguments:
% 
%                            im: Input Image
%         
%                         sigma: standard deviation of the gaussian
% 
%                                Caution: If the pixel spacing is not 
%                                specified then unit isotropic spacing is 
%                                assumed as a result of which the units of 
%                                sigma will be assumed to be pixels.
% 
%                     direction: Specifies the dimension along which to
%                                compute the gradient.
% 
%                                Must be an integer between 1 and ndims(im)
%                               
%                                Note: direction = 1 gives the gradient 
%                                along the y-axis and direction = 2 gives 
%                                the gradient along the x-axis.
% 
%     Optional Input Arguments:
% 
%                       spacing: pixel spacing of the input image.
%                                can be a scalar value or an array of size
%                                equal to the number of image dimensions.
%                                Default: 1 (isotropic spacing)                      
%                 
%               borderCondition: specifies the way in which the image is 
%                                padded at the borders. This argument is
%                                supplied as input to the functuin 'padarrayXT'. 
%                                Default: 'symmetric'
%                                Options: 'symmetric', 'replicate', 'circular', 
%                                         'antisymmetric', or a constant value
% 
%      UseNormalizedDerivatives: true/false
%                                specifies whether the gaussian derivatives
%                                should be scale-normalized.
%                                Default: false
% 
%         UseNormalizedGaussian: true/false
%                                specifies whether the guassian kernel should
%                                be normalized or not
%                                Default: true
%                                   
%                        UseGPU: true/false
%                                A flag that specifies whether or not to
%                                use the GPU for convolution.  
%                                   True - Uses GPU if installed
%                                   False - otherwise (default-value)
%                                   
%     Output Arguments:
% 
%                      imGaussGrad: gaussian filtered image
%         
%     Example-1: Application of Derivative of Gaussian Filter to a 1D Step Edge
%         
%         step = @(x) ( double( x >= 0 ) ); % 1-D step edge
%         g = @(x,sigma) ( (sigma * sqrt(2*pi))^-1 * exp( -x.^2 / (2*sigma^2) ) ); % gaussian
%         truedogresponse = @(x,sigma) ( g(x,sigma) ); % true LoG response
%         
%         sampleSpacing = 0.01;
%         x = -20:sampleSpacing:20;
%         
%         % check the difference between the response of this filter and the true LoG response
%         figure; hold all;
%         sigma = 2;
% 
%             % plot step edge
%             plot( x, step(x), '-' );
% 
%             % plot true scale-normalized (multiply by sigma) DOG response
%             plot( x, sigma * truedogresponse(x,sigma), '-' ); 
% 
%             % plot response of discrete DOG response 
%             plot( x, filterGaussGradND( step(x), sigma, 1, ...
%                                         'spacing', sampleSpacing, ...
%                                         'UseNormalizedDerivatives', true ), '-' ); 
% 
%             title( 'Comparison of our filter with true DOG response for a 1D Sted Edge', ...
%                    'FontWeight', 'bold' );
%             legend( '1D Step Edge', ... 
%                      sprintf( 'True DOG Response (\\sigma = %d)', sigma ), ...
%                     'Our  DOG Response' );          
% 
%         % check the difference between normalized and unnormalized gaussian derivatives
%         figure;
%         sigmaValues = 2.^[0:2];
% 
%             % plot responses for unnormalize gaussian derivatives
%             subplot(2,1,1); hold all;
%             strLegend = {};
%             
%                 % plot step edge - scale to reduce height of ridge so the tiny
%                 % unnormalized DOG responses become visible
%                 plot( x, 0.01 * step(x), '-' ); 
% 
%                 strLegend{end+1} = sprintf( '0.01 \\times 1D Step Edge' );
%                 
%                 % plot response of discrete DOG response 
%                 for i = 1:numel( sigmaValues )
%                     
%                     plot( x, filterGaussGradND( step(x), sigmaValues(i), 1, ...
%                                                 'spacing', sampleSpacing, ...
%                                                 'UseNormalizedDerivatives', false ), '-' ); 
%                                        
%                     strLegend{end+1} = sprintf( 'DOG Response for \\sigma = %d', sigmaValues(i) );
% 
%                 end
% 
%             title( 'DOG Responses of a 1D Step Edge for Unnormalized Gaussian Derivatives', ...
%                    'FontWeight', 'bold' );
%             legend( strLegend );          
% 
%             % plot responses for scale-normalized gaussian derivatives
%             subplot(2,1,2); hold all;
%             strLegend = {};
%             
%                 % plot step edge
%                 plot( x, step(x), '-' );
% 
%                 strLegend{end+1} = sprintf( '1D Step Edge' );
%                 
%                 % plot response of discrete DOG response 
%                 for i = 1:numel( sigmaValues )
%                     
%                     plot( x, filterGaussGradND( step(x), sigmaValues(i), 1, ...
%                                                 'spacing', sampleSpacing, ...
%                                                 'UseNormalizedDerivatives', true ), '-' ); 
%                                        
%                     strLegend{end+1} = sprintf( 'DOG Response for \\sigma = %d', sigmaValues(i) );
% 
%                 end
% 
%             title( 'DOG Responses of a 1D Step Edge for Normalized Gaussian Derivatives', ...
%                    'FontWeight', 'bold' );
%             legend( strLegend );          
%  
%     Author: Deepak Roy Chittajallu
% 
% 

    p = inputParser;
    p.CaseSensitive( false );
    p.addRequired( 'im', @(x) ( isnumeric(x) ) ); 
    p.parse( im );
    
    imdims = ndims(im);   
    if numel(im)==max(size(im))
        imdims = 1;
    end           
    
    p.addRequired( 'sigma', @(x) ( isnumeric(x) && (isscalar(x) || numel(sigma) == imdims) ) ); 
    p.addRequired( 'direction', @(x) (isscalar(x) && ismember(x, 1:imdims)) ); 
    p.addParamValue( 'spacing', ones( 1, ndims(im) ), @(x) (isnumeric(x) && ismember(numel(x), [1,imdims])) );
    p.addParamValue( 'borderCondition', 'symmetric', @(x) ( isscalar(x) || (ischar(x) && ismember( lower(x), { 'symmetric', 'replicate', 'antisymmetric' } ) ) ) );
    p.addParamValue( 'UseNormalizedDerivatives', false, @(x) (isscalar(x) && islogical(x)) );
    p.addParamValue( 'UseNormalizedGaussian', true, @(x) (isscalar(x) && islogical(x)) );
    p.addParamValue( 'UseGPU', false, @(x) (isscalar(x) && islogical(x)) );
    p.parse( im, sigma, direction, varargin{:} );
    
    spacing = p.Results.spacing;    
    borderCondition = p.Results.borderCondition;
    flagNormalizeDerivatives = p.Results.UseNormalizedDerivatives;
    flagNormalizeGaussian = p.Results.UseNormalizedGaussian;
    flagUseGPU = p.Results.UseGPU;
    
    % adjust sigma according to spacing
    sigma = sigma / spacing(direction);
    
    % Compute Derivative of Gaussian Kernel
    w = ceil(3 * sigma);
    x = -w:w;
    
    gaussKernel = exp(-x.^2 / (2*sigma^2));    
    if flagNormalizeGaussian
        gaussKernel = gaussKernel / sum(gaussKernel(:));
    end
    
    dogKernel = -(x / sigma^2) .* gaussKernel;
    if flagNormalizeDerivatives
        dogKernel = sigma * dogKernel;
    else
        
    end
    dogKernel = dogKernel - mean( dogKernel(:) ); % must sum to zero
    
    if imdims > 1
        kernelShape = ones(1,imdims);
        kernelShape(direction) = numel(dogKernel);    
        dogKernel = reshape( dogKernel, kernelShape );
    end
    
    % apply kernel by convolution
    if imdims > 1
        padsize = zeros(1,imdims);
        padsize(direction) = w;
    else
        padsize = [0, w];
    end
    
    imPadded = padarrayXT(im, padsize, borderCondition);
    imGaussGrad = convnfft( imPadded, dogKernel, 'valid', 'UseGPU', flagUseGPU);
    
end