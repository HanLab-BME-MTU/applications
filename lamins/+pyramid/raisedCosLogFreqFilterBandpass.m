function [ bandpass ] = raisedCosLogFreqFilterBandpass( r)
%raisedCosLogFreqFilter implements a raised cosine filter in the frequency
%domain

lowerLimit = pi/4;
upperLimit = pi;

lowMask  = abs(r) <= lowerLimit;
highMask = abs(r) >= upperLimit;

transitionMask = ~highMask & ~lowMask;
logr = log2(r(transitionMask) * 2/pi);

bandpass = zeros(size(r));
bandpass(transitionMask) = cos(pi/2*logr);


end

