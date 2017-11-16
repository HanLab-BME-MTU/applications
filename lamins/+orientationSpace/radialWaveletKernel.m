function [ radialFilter ] = radialWaveletKernel( f_c, b_f, N )
%radialWaveletKernel Calculate the orientationSpace radial wavelet kernel
%
% See also orientationSpace.kernel, orientationSpace.angularKernel,
% orientationSpace.getFrequencySpaceCoordinates

% Mark Kittisopikul
% Jaqaman Lab
% UT Southwestern
% June 24th, 2016

if(nargin < 2 || isempty(b_f))
    b_f = f_c/sqrt(2);
end
if(nargin < 3)
    N = 1024;
end

assert(length(f_c) == length(b_f),'orientationSpace.radialKernel:f_c and b_f must be of equal size');

f_c = shiftdim(f_c(:),-3);
b_f = shiftdim(b_f(:),-3);

if(isstruct(N))
    coords = N;
else
    coords = orientationSpace.getFrequencySpaceCoordinates(N);
end

%% Radial part
% compute radial order, f_c = sqrt(K_f * b_f^2)
if(f_c)
    K_f = (f_c ./ b_f).^2;
% 
%     % scale frequency
%     % f_s = coords.f / f_c;
%     f_s = bsxfun(@rdivide,coords.f,f_c);
% 
%     % Equation 3.11
%     % Note -(f^2 - f_c^2)/(2*b_f^2) = (1 - (f/f_c)^2)/(2* b_f^2/f_c^2)
%     %                               = (1 - (f/f_c)^2)/(2 / K_f)
%     % radialFilter = f_s.^K_f .* exp((1 - f_s.^2)*K_f/2);
%     radialFilter = bsxfun(@power,f_s,K_f);
%     radialFilter = radialFilter .* exp(bsxfun(@times,(1 - f_s.^2),K_f/2));
%     % radialFilter2 = f_s^K_f .* exp(-(f.^2-f_c.^2)/2/b_f.^2);
%     % assertEqual(radialFilter,radialFilter2);


    r = coords.f / f_c * pi/2;
    
    %% From pyramid.raisedCosLogFreqFilterBandpass where maximum occurs at
    % pi/2. Limit are at pi/4 to pi
    lowerLimit = pi/8*2^sqrt(K_f);
    upperLimit = 2*pi*2^-sqrt(K_f);

    lowMask  = abs(r) <= lowerLimit;
    highMask = abs(r) >= upperLimit;

    transitionMask = ~highMask & ~lowMask;
    % logr ranges from -1 to 1 so that pi/2*logr ranges from -pi/2 to pi/2
    % max will occur when logr == 0 (and thus when r = pi)
    logr = log2(r(transitionMask) * 2/pi);

    bandpass = zeros(size(r));
    bandpass(transitionMask) = cos(pi/2*logr/(2-sqrt(K_f)));
    
    
    %% End excerpt
    radialFilter = bandpass;
else
%     radialFilter = exp(-bsxfun(@rdivide,coords.f.^2,2*b_f.^2));
%     radialFilter(1) = 0;
    
    % This should be a low pass filter based on b_f
end

end

