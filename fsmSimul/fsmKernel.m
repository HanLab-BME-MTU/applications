function [M,oPoly,pixPicPolym] = fsmKernel(mic,iPoly,nRepl)
%FSMKERNEL kernel of the FSM (fluorescent speckle microscopy) simulation  
%
% SYNOPSIS [M,oPoly,pixPicPolym] = fsmKernel(mic,iPoly,nRepl)
%
% INPUT  mic parameter structure of the microscope containing the fields
%          .monomSze : monomer size
%          .LR : labeling ratio
%          .NA : numerical aperture 
%          .mag: magnification
%          .wvl: exitation wavelength
%          .pixSize: pixel size in [nm]
%          .polymSize : length of the polymer in [nm]
%          .monomSize : length of a monomer in [nm]
%        
%        To get a full version of the data structure (even more than what is currently
%        needed by fsmKernel call 
%                      mic = getFsmDflts();
%
%        iPoly : input polymer (if none exists set ipoly = []; in this case nRepl is 
%                is ignored and a polymer with nOfPolyms = .polymSize / .monomSize is generated)
%        nRepl : number of monomers to be added at the front of iPoly and 
%                to be removed at the end (assuming a treadmilling, polymer with constant 
%                length) 
%
% OUTPUT M     : modulation as a measure of the quality of the speckle signal
%        oPoly : output polymer (iPoly but with nRepl replaced monomers)
%        pixPicPolym : pixelated picture of the polymer
%
% SEE ALSO getFsmDflts()   for more information on the data structure
%        
% started Nov-15, 2000 : gD

% some data security checking
if isempty(iPoly)
    iPoly = zeros(round(mic.polymSize/mic.monomSize),1);
    nRepl = round(mic.polymSize/mic.monomSize);
end;
if(length(iPoly)<nRepl)
   error('number of monomers to be replaced larger than polymer');
end;

% calculate how many monomers fall inside one pixel
nPix = round(mic.pixSize / (mic.mag * mic.monomSize));
                              % how many monomers do fall inside a pixel

% adding monomers at the front, cutting the same number at the end
% Labelled monomers have a 1, unlabeled a 0
oPoly = [(rand(nRepl,1)<mic.LR);iPoly(1:end-nRepl)];

% add some some stupid, additive noise 
addNPoly = 1/100 * rand(length(oPoly),1);
multNPoly = 1/2 * rand(length(oPoly),1);
totPoly = oPoly + oPoly.*multNPoly + addNPoly;

% convolve with PSF and sample in pixels
	   
% 1.Step: create PSF kernel
%  defintion between -/+ the 2nd root of the Bessel function 
y = -7.02*(1.22 * mic.wvl/(2*mic.NA))/3.83 + 0.01 : mic.monomSize : ...
   7.02*(1.22 * mic.wvl/(2*mic.NA))/3.83;
ys = y /(1.22 * mic.wvl/(2*mic.NA))*3.83;
psfs = (besselj(1,ys)./ys);
psf = psfs.*conj(psfs);
psf = psf / sum(psf);  % make sure that the optics is an energy 
                       % preserving filter

% 2. Step: convolution
picPolym = conv(totPoly,psf);

% generate pixelated picture
pixBdIndx = 1:nPix:length(picPolym);
cumPicPolym = cumsum(picPolym);
% clear pixPicPolym;
for(i=2:length(pixBdIndx))
   pixPicPolym(i-1)= cumPicPolym(pixBdIndx(i))- ...
      cumPicPolym(pixBdIndx(i-1));
end;

% compute a contrast measure
M = fsmGetContrastMeas('meanModulation',pixPicPolym);

