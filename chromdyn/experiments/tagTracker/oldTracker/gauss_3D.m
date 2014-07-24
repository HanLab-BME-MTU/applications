function gauss=gauss_3D(amp,sigma,cent,XX,YY,ZZ)
% wrapper to ensure backwards compatibility. gauss_3D will be discontinued. Use gaussListND instead.

gauss = amp*GaussListND([XX,YY,ZZ],sigma,cent,0);