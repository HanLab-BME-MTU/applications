function [sizeXY, sizeZ] = calcFilterParms(wvl, NA, n, gaussOrBessel, corr, pixelSize)
% calcFilterParms calculates the width of the psf or the gaussian approximation from microscope parameters
%
% SYNOPSIS filterparms = ...
%               calcFilterParms(wvl, NA, n, gaussOrBessel, corr,pixelsize)
%
% INPUT    wvl:  (opt) wavelength in microns {0.525}
%          NA :  (opt) numerical aperture {1.4}
%          n  :  (opt) refractive index of oil/air {1.51}
%          gaussOrBessel : (opt) whether to calculate gaussian or bessel
%                           psf. {'gauss'} / 'bessel'
%          corr : (opt) [corrXY, corrZ]. Correction factors. {1,1}
%          pixelSize : (opt) [sizeX, sizeZ]. If size is given, the
%                            psf-width is returned in pixel
%
% OUTPUT   siyeXY, sizeZ: Size of psf in microns or pixels
%
% c: jonas, 10/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=======================
% DEFAULTS & TEST INPUT
%=======================

def_wvl = 0.525;
def_NA  = 1.4;
def_n   = 1.518;
def_gob = [0.21, 0.66];
def_corr= [1,1];
def_pix = [1,1];

% not checking for weird input numbers.
if nargin < 1 || isempty(wvl)
    wvl = def_wvl;
end
if nargin < 2 || isempty(NA)
    NA = def_NA;
end
if nargin < 3 || isempty(n)
    n = def_n;
end
if nargin < 4 || isempty(gaussOrBessel)
    gob = 'def_gob';
else
    switch lower(gaussOrBessel)
        case 'gauss'
            % multiplication factor for Gauss fitted to Bessel (see Dom's
            % first JM paper)
            gob = [0.21, 0.66];
        case 'bessel'
            % multiplication factor for Bessel (see any microscopy book)
            gob = [0.61, 2];
        otherwise
            error('non recognized option for ''gaussOrBessel''')
    end
end
if nargin < 5 || isempty(corr)
    corr = def_corr;
end
if nargin < 6 || isempty(pixelSize)
    pixelSize = def_pix;
end

%=======================


%=======================
% CALCULATE PSF-SIZE
%=======================
    
sizeXY = corr(1)*(gob(1)*wvl/NA)/pixelSize(1);
sizeZ  = corr(2)*(gob(2)*n*wvl/NA^2)/pixelSize(2);
