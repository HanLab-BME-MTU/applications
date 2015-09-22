function [ imSmooth ] = filterGaussND(im, sigma, varargin )
% filterGaussND: An ND implementation of gaussian smoothing
% 
%     [ imSmooth ] = filterGaussND(im, sigma, varargin );
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
%                                   
%     Output Arguments:
% 
%                      imSmooth: gaussian filtered image
%         
%     Author: Deepak Roy Chittajallu
% 

    p = inputParser;
    p.CaseSensitive( false );
    p.addRequired( 'im', @(x) ( isnumeric(x) ) ); 
    p.parse( im );
    
    p.addRequired( 'sigma', @(x) ( isnumeric(x) && (isscalar(x) || numel(sigma) == ndims(im)) ) ); 
    p.addParamValue( 'spacing', ones( 1, ndims(im) ), @(x) (isnumeric(x) && ismember(numel(x), [1,ndims(im)])) );
    p.addParamValue( 'borderCondition', 'symmetric', @(x) ( isscalar(x) || (ischar(x) && ismember( lower(x), { 'symmetric', 'replicate', 'antisymmetric' } ) ) ) );
    p.addParamValue( 'UseNormalizedGaussian', true, @(x) (isscalar(x) && islogical(x)) );
    p.addParamValue( 'UseGPU', false, @(x) (isscalar(x) && islogical(x)) );
    p.parse( im, sigma, varargin{:} );
    
    spacing = p.Results.spacing;    
    borderCondition = p.Results.borderCondition;
    flagNormalizeGaussian = p.Results.UseNormalizedGaussian;
    flagUseGPU = p.Results.UseGPU;
    
    imdims = ndims( im );  
    
    if isscalar( sigma )
        sigma = sigma + zeros(1,imdims);
    end
    
    % adjust sigma according to spacing
    sigma = sigma ./ spacing;
    
    % Compute Gaussian Kernel
    w = ceil(3 * sigma);
    xrange = cell(1,imdims);
    
    for i = 1:imdims
       xrange{i} = -w(i):w(i); 
    end
    
    x = cell(1,imdims);
    if imdims > 1
        [x{:}] = ndgrid( xrange{:} );    
    else
        x{1} = xrange{1};
    end
        
    gaussKernel = ones( size(x{1}) );
    for i = 1:imdims
        gaussKernel = gaussKernel .* ( exp(-x{i}.^2 / (2*(sigma(i))^2)) );
    end        
    
    if flagNormalizeGaussian
        gaussKernel = gaussKernel / sum(gaussKernel(:));
    end
    
    % apply guassian kernel
    if imdims > 1
        padsize = w;
    else
        padsize = [0, w];
    end    
    imPadded = padarrayXT(im, padsize, borderCondition);
    imSmooth = convnfft( imPadded, gaussKernel, 'valid', 'UseGPU', flagUseGPU);
    
end
