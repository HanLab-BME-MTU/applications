function [rayleighLimit, angle] = rayleighFromOri(normedVectors, wavelength, NA, n)
%RAYLEIGHFROMORI calculates the Rayleigh limit given a vector in space
%
% SYNOPSIS: [rayleighLimit, rayleighRatio, angle] = rayleighFromOri(normedVectors, dataProperties)
%
% INPUT normedVectors: n-by-3 list of normed vectors
%       wavelength (opt): wavelength in microns {0.525}
%       NA (opt): numerical aperture {1.4}
%       n (opt): refractive index of immersion oil {1.518}
%
% OUTPUT rayleighLimit: rayleighLimit in the specified orientation (in
%           microns)
%	     angle: angle of the vector with the xy-plane
%
% REMARKS to calculate the ratio of the vectorLength to the rayleighLimit,
%           simply do vectorLength./rayleighLimit
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 09-May-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%========================
%% TEST INPUT
%========================

% defaults
def_wvl = 0.525;
def_NA  = 1.4;
def_n   = 1.518;

if nargin < 1 || isempty(normedVectors)
    % don't error, as this is a helper function
    rayleighLimit = [];
    angle = [];
    return
end

% check for optional input.
if nargin < 2 || isempty(wavelength)
    wavelength = def_wvl;
end
if nargin < 3 || isempty(NA)
    NA = def_NA;
end
if nargin < 4 || isempty(n)
    n = def_n;
end

%========================


%===========================
%% CALCULATE RAYLEIGH-LIMIT
%===========================

% orientation: [-pi/2...pi/2]. angle = acos(v*e_z/(|v|*|e_z|)). Because the
% vector is normed, all that remains in the brackets is the third component
% of the normedVectors

% norm the vectors
angle = pi/2-acos(normedVectors(:,3));

% calculate Rayleigh limit

[ralXY, ralZ] = calcFilterParms(wavelength, NA, n, 'bessel');

rayleighLimit = sqrt(1./ (cos(angle).^2 / ralXY^2 + ...
    sin(angle).^2 / ralZ^2) );

