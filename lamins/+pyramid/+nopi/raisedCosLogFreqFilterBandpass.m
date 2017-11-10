function [ bandpass ] = raisedCosLogFreqFilterBandpass( r)
%raisedCosLogFreqFilter implements a raised cosine filter in the frequency
%domain

% no pi so [-0.5,0.5] range corresponds to be [-pi,pi]

% achieves absolute maximum value of 1 at 0.25, corresponding to pi/2

lowerLimit = 1/8; % pi/4 / (2*pi)
upperLimit = 1/2; % pi / (2*pi)

lowMask  = abs(r) <= lowerLimit;
highMask = abs(r) >= upperLimit;

transitionMask = ~highMask & ~lowMask;
logr = log2(r(transitionMask)*4);

bandpass = zeros(size(r));
bandpass(transitionMask) = cos(pi/2*logr);


end

