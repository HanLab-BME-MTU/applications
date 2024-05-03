function [] = erodeForceField( MD,band )
%function [] = erodeForceField( MD,band ) imports the processed forceField
%from step 4 of TFM package erode the boundary force coeeficient
% input:    pathForTheMovieDataFile:    path to the movieData file or MD
%                                       itself
%           band:                       band width for cutting border
%                                       (default=4)
%           tmax:                       maximum traction, if not indicated,
%                                       80% of the true global maximum will
%                                       be used.
% 
% output:   
%           strainEnergy: 1/2*integral(traction*u) in femtoJ (10^-15 Joule)
%           images of heatmap stored in pathForTheMovieDataFile/heatmap
if nargin <2
    band = 1;
end
