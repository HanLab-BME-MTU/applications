function [ maxima_splines, minima_splines, response ] = splinesForExtremaFT( response, K, freq )
%splinesForExtremaFT Obtain spline fits for the local maxima through K

if(nargin < 3)
    freq = false;
end

%% Expand response
if(isvector(response))
    response = response(:);
    response = [response getResponseAtOrderFTatPoint(response(:),K(1),K(2:end),freq)];
else
end

[maxima,minima] = interpft_extrema(response);

maxima_splines = getSplineFromExtrema(maxima,response,K);
if(nargout > 1)
    minima_splines = getSplineFromExtrema(minima,response,K);
end

end

function response = getResponseAtOrderFTatPoint(r,K_old,K_new,freq)
    % Get response at a lower order using Fourier Transform at a particular point
    % INPUT
    % r - row, scalar integer
    % c - column, scalar integer
    % K_new - new angular order
    % OUTPUT
    % response at pixel (r,c) at angular order K_new

    % New number of coefficients
%             n_new = 2*ceil(K_new)+1;
    n_old = 2*ceil(K_old)+1;
    n_new = 2*K_new+1;

    % The convolution of two Gaussians results in a Gaussian
    % The multiplication of two Gaussians results in a Gaussian
    % The signal has been convolved with a Gaussian with sigma = pi/obj.n
    % Compute the signal convoled with a Gaussian with sigma = pi/n_new
    % Note if n == n_new, we divide by 0. Then s_inv = Inf
    s_inv = sqrt(n_old^2.*n_new.^2./(n_old.^2-n_new.^2));
    s_hat = s_inv/(2*pi);
    x = -ceil(K_old):ceil(K_old);

    % Each column represents a Gaussian with sigma set to s_hat(column)
    f_hat = exp(-0.5 * bsxfun(@rdivide,x(:),s_hat).^2); % * obj.n/n_new;
    f_hat = ifftshift(f_hat,1);

    % Angular response will be in a column
    if(freq)
        a_hat = squeeze(r);
    else
        a_hat = fft(squeeze(real(r)));
    end
    % Each column represents an angular response with order K_new(column)
    a_hat = bsxfun(@times,a_hat,f_hat);
    response = ifft(a_hat);
end

function sp = getSplineFromExtrema(extrema,response,K)
    nDerivs = 5;
    extrema = orientationSpace.diffusion.alignExtrema(extrema);
    extrema = fliplr(extrema);
    K = fliplr(K);
    response = fliplr(response);
    validTracks = sum(~isnan(extrema),2) > 1;
    extrema = extrema(validTracks,:);
    derivs = orientationSpace.diffusion.orientationMaximaDerivatives(response,K,nDerivs,extrema);
    
    

    
    x = joinColumns(repmat(K,nDerivs+1,1));
    for trackNum = 1:size(extrema,1)
        y = joinColumns([extrema(trackNum,:); permute(derivs(trackNum,:,:),[3 2 1])]);
        valid = ~isnan(y);
        knots = optknt(x(valid).',nDerivs+3);
        sp(trackNum) = spapi(knots,x(valid).',y(valid).');
    end
end