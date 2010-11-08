% runLocalization performs Gaussian-approximated PSF fitting for all data sets in the input structure.
%
% INPUT     data        : array of experiment structures
%           sigma       : standard deviation of the Gaussian PSF approximation.
%                         Compute using 'getGaussianPSFsigma.m'
%           {overwrite} : 

% Francois Aguet, October 2010

function runLocalization(data, overwrite)

if nargin<2
    overwrite = 0;
end

nExp = length(data);
parfor i = 1:nExp   
    psfLocalization(data(i), overwrite);
end