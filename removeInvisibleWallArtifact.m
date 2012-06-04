function mask = removeInvisibleWallArtifact(mask,binning)
%REMOVEINVISIBLEWALLARTIFACT removes the artifact in some of bobs confocal images 
% 
% mask = removeInvisibleWallArtifact(mask,binning)
% 
% This removes the cylindrical artifact created in the segmentation by a
% problem with one of the pinholes in bob's CSUX confocal head. It's a
% quick hack based on some blank correction images he took.
% 
% 
% 
% Hunter Elliott
% 5/8/2012
%


[M,N,P] = size(mask);

showPlots = false;%Plots for development/debugging

%NOTE: I guess this could be replaced with standard correction masks for
%each image size, but given orchestra transfer speeds generating the masks
%on the fly is probably faster...


%Parameters & function for quadratic which was fit to artifact locations
%extracted from correction images. Technically this should be a circle, but
%this gave a more-than-good-enough fit given it is a small arc of a circle.

quadPar = [0.000286858841885258,0.0108763662090231,73.3785278830350];
quadFun = @(x,p)(p(1)*x.^2 + p(2)*x + p(3));

%Find pixels in the image which match this quadratic
[X,Y] = meshgrid(-N/2+1:N/2,M:-1:1);
X = X .* (binning/2);
Y = Y .* (binning/2);
Xquad = quadFun(X,quadPar);
%Rounding both ways to give a little leeway to cover the entire artifact,
%and repeat over z
corrMask = repmat(ceil(Xquad) == round(Y) | floor(Xquad) == round(Y),[1 1 P]);

if showPlots
    figure
    x = -N/2:N/2;
    y = M - quadFun(x,quadPar);
    %imshow(sum(mask,3),[])
    imshow(mean(double(image),3),[])
    hold on
    plot(x+N/2,y,'-')    
    spy(corrMask(:,:,1),'r')
        
end

%Now remove these voxels from the mask. This will create a 1-2 pixel slice
%through any part of the cell which passes through the wall, but the
%post-processing will handle this just as it does holes due to noise etc.
mask(corrMask(:)) = false;




