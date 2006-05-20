function [specDenV,lambda,errFlag] = spectralDensity(gamma)
%SPECTRALDENSITY calculates the spectral density from the autocovariance function
%
%SYNOPSIS [specDenV,lambda,errFlag] = spectralDensity(gamma)
%
%INPUT  gamma   : Autocovariance function of series or process, 
%                 where gamma(i) is the autocovariance at lag i-1.
%
%OUTPUT specDenV: Spectral density.
%       lambda  : Values of lambda for calculating the spectral density.
%       errFlag : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, March 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spectralDensityV = [];
lambda = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if correct number of arguments were used when function was called
if nargin ~= nargin('spectralDensity')
    disp('--spectralDensity: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Spectral Density Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%determine values of lambda for spectral density calculation
lambda = [0:0.02:pi]';
numLambda = length(lambda);

%reserve memory for specDenV
specDenV = zeros(numLambda,1);

%find the maximum lag whose autocorrelation is supplied
maxLag = length(gamma) - 1;

%calculate the spectral density
for i=1:numLambda
    specDenV(i) = (2*cos([1:maxLag]*lambda(i))*gamma(2:end) ...
        + gamma(1))/2/pi;
end


%%%%% ~~ the end ~~ %%%%%
