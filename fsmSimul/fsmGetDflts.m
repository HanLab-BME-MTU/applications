function data = fsmGetDflts
%FSMGETDFLTS gets a reasonable set of parameters controling FSM
%  
% SYNOPSIS   data = fsmGetDflts
% 
% INPUT      none
%
% OUTPUT     data is a structure with the following parameters
%            *.polymSize = 10000;  (size of the polymer)
%            *.monomSize = 5;      (size of the monomer)
%            *.pixSize = 12000;    (pixel size in the image plane)
%            *.mag = 100;          (microscope magnification)
%            *.NA  = 1.4;          (numerical aperture)
%            *.wvl = 500;          (wavelength)
%            *.LR = 5/100;         (percentage of labelled free monomers)
%            *.camRng = [0,4095];  camera range in greyvalues



% definition of key quantities (units nm)
data.polymSize = 50000;  % (size of the polymer)
data.monomSize = 5;      % (size of the monomer)
data.pixSize = 6700;     % (pixel size in the image plane)
data.mag = 100;          % (microscope magnification)
data.NA  = 1.4;          % (numerical aperture)
data.wvl = 500;          % (wavelength)
data.LR = 0.5/100;         % (percentage of labelled free monomers)
data.camRng = [0,4095];  % camera range in greyvalues