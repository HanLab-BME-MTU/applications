function[avmod]=averageSpeckleModulation(cands)
%calculates average speckle modulation in image from cands
%SYNOPSIS [avmod]=averageSpeckleModulation(cands)
%speckles = positions of significant maxima = maxima with status==1
%input : cands
%output: avmod = average modulation of all speckles in cand
%
mstat=cat(1,cands.status);
disp(['no. of speckles',num2str(max(size(nonzeros(mstat))))]);
%maximum intensity vector (for all maxima in plane)
imax=cat(1,cands.ILmax);
%average intensity of three surr. minima (for all maxima in plane)
imin=cat(1,cands.IBkg);
%modulation of counted speckles
mod=nonzeros(mstat.*(imax-imin)./(imax+imin));
%average
avmod=mean(mod);