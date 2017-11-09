function [ lowpass, highpass, bandpass ] = raisedCosLogFreqFilter( r)
%raisedCosLogFreqFilter implements a raised cosine filter in the frequency
%domain

lowerLimit = pi/4;
upperLimit = pi/2;

lowMask  = r <= lowerLimit;
highMask = r >= upperLimit;

transitionMask = ~highMask & ~lowMask;
logr = log2(r(transitionMask) * 2/pi);

lowpass = double(lowMask);
% lowpass(lowMask) = 1;
lowpass(transitionMask) = -sin(pi/2 * logr);
% lowpass(highMask) = 0;

if(nargout > 1)
    highpass = double(highMask);
%     highpass(lowMask) = 0;
    highpass(transitionMask) = cos(pi/2 * logr);
%     highpass(highmask) = 1;
end

if(nargout > 2)
    bandpass = lowpass.*highpass;
end


end

