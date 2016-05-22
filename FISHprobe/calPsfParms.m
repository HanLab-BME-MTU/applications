function [psfSigma, filterParms] = calPsfParms(chaNum, dataProperties)
% calPsfParms Prepare parameters for psf filter and image deconvolution
% channel
%
% INPUT:
%   chaNum: 
%   dataProperties
%
% OUTPUT:
%   psfSigma:
%   filterParms:
%
% 05/2016 Ning Zhang


% Excitation or emission???
wvl = dataProperties.channel(chaNum).emissionWavelength;

p = inputParser;
p.addRequired('wvl', @(x) (isnumeric(x) && ~isempty(x)));
p.addRequired('NA', @(x) (isnumeric(x) && ~isempty(x)));
p.addRequired('n', @(x) (isnumeric(x) && ~isempty(x)));
p.addRequired('pixelSize', @(x) (isnumeric(x) && ~isempty(x)));
% sigmaCorrection defined by default
p.addOptional('corr', [1,1], @(x) (isnumeric(x) && ~isempty(x)));
p.addParameter('gaussOrBessel', 'gauss', @isstr);
p.addParameter('detectionMethod', 'mnp', @isstr);

p.parse(wvl, dataProperties.NA, dataProperties.refractiveIndex, ...
    [dataProperties.PIXELSIZE_XY, dataProperties.PIXELSIZE_Z]);

wvl = p.Results.wvl;
NA = p.Results.NA;
n = p.Results.n;
pixelSize = p.Results.pixelSize;
corr = p.Results.corr;
gaussOrBessel = p.Results.gaussOrBessel;
switch lower(gaussOrBessel)
    case 'gauss'
        % multiplication factor for Gauss fitted to Bessel (see Dom's first JM paper)
        gob = [0.21, 0.66];
    case 'bessel'
        % multiplication factor for Bessel (see any microscopy book)
        gob = [0.61, 2];
    otherwise
        error('non recognized option for ''gaussOrBessel''')
end

% Calculate psf sigma on pixel level (Refer to calcFilterParms.m by Jonas)
sigmaXY = corr(1)*(gob(1)*wvl/NA)/pixelSize(1);
sigmaZ  = corr(2)*(gob(2)*n*wvl/NA^2)/pixelSize(2);
% Calculate patch size on pixel level (Refer to defaultDataProperties.m by Jonas)
patchXYZ = roundOddOrEven(4*[sigmaXY sigmaXY sigmaZ], 'odd', 'inf');
filterParms = [sigmaXY, sigmaXY, sigmaZ, patchXYZ];
psfSigma = [sigmaXY, sigmaXY, sigmaZ];
