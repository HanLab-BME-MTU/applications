function simIm = generateSynthetic3DCellImage(varargin)



showPlots = true;

%% -------------- Input -------------- %%

ip = inputParser;
ip.FunctionName = mfilename;

%Image parameters
ip.addParamValue('ImageSize',[520 696 60],@(x)(numel(x)==3 && all(x>1)));
ip.addParamValue('ZSpacing',500e-9,@(x)(numel(x)==1 && x > 0));
ip.addParamValue('ImageType','Sphere',@(x)(ischar(x)));
ip.addParamValue('MembraneIntensity',250,@(x)(numel(x)==1 && x > 0));

%Background Parameters
ip.addParamValue('NoiseLambda',50,@(x)(x>0 && numel(x) == 1));
ip.addParamValue('BackgroundOffset',150,@(x)(x>0 && numel(x) == 1));
%PSF Parameters, input
ip.addParamValue('DistFromCoverslip',200e-6,@(x)(x>0 && numel(x) == 1));%Distance of top image z-plane from sample side of coverslip. Bob estimates 25-250um in practice.
ip.addParamValue('Magnification',60,@(x)(x>0 && numel(x) == 1));
ip.addParamValue('PsfKernelSize',12);

ip.parse(varargin{:});
p = ip.Results;

%% -------------- Fixed Parameters -------------- %%
%Designed to replicate images from Bob's spinning disk


%PSF Parameters
psfProp.Ti0= 2.7000e-4;
psfProp.Ni0= 1.3300;
psfProp.Ni= 1.3300;
psfProp.Tg0= 1.7000e-4;
psfProp.Tg= 1.7000e-4;
psfProp.Ng0= 1.5150;
psfProp.Ng= 1.5150;
psfProp.Ns= 1.4300;
psfProp.lambda= 581e-9;%Lambda max for tdTomato FP
psfProp.M= p.Magnification;
psfProp.NA= 1.2000;
psfProp.alpha= 1.1250;
psfProp.pixelSize= 6.4000e-06;

%Background parameters
noiseLambda = 50;%Lambda of background poisson noise
backOffset = 150; %Camera intensity offset.





%% --------------- Image Simulation -------------- %%


   % ------------- Image Geometry Generation --------- %
   
   
switch p.ImageType
    
    case 'Sphere'
        
        [X,Y,Z] = meshgrid(-p.ImageSize(2)/2:p.ImageSize(2)/2-1,-p.ImageSize(1)/2:p.ImageSize(1)/2-1,-p.ImageSize(3)/2:p.ImageSize(3)/2-1);
        R = sqrt( X .^2 + Y .^2 + Z .^2);
        imageGeo = zeros(p.ImageSize);
        imageGeo(R < 20 & R > 18) = p.MembraneIntensity;
    
end



   % ------------- PSF Generation & Convolution --------- %

%Generate the simulated PSF
%TEMP - should use different x,y,z for different image locations?? How much
%does it matter?
zPlanePos = p.DistFromCoverslip - ((p.PsfKernelSize-1) * p.ZSpacing) : p.ZSpacing : p.DistFromCoverslip + ((p.PsfKernelSize-1) * p.ZSpacing); 

psfEmit = vectorialPSF(0,0,p.DistFromCoverslip,zPlanePos,p.PsfKernelSize,psfProp);
%Now get the excitation PSF
psfProp.lambda = 560e-9;%Laser line bob uses for illumination. Lambda max excitation is 554nm.
psfExci = vectorialPSF(0,0,p.DistFromCoverslip,zPlanePos,p.PsfKernelSize,psfProp);

%Combine emission and excitation PSFs to get effective confocal PSF
psfEff = psfEmit .* psfExci;
psfEff = psfEff ./ sum(psfEff(:));%Normalize for use in convolution.

if showPlots
    psfFig = figure;    
    subplot(1,2,1);
    imagesc(squeeze(max(psfEff,[],1))),axis image,colormap gray
    title('PSF, yz')
    subplot(1,2,2);
    imagesc(squeeze(max(psfEff,[],3))),axis image,colormap gray    
    title('PSF, xy')
end

%Convolve PSF with input image
simIm = convn(padarrayXT(imageGeo, p.PsfKernelSize * ones(1,3),'symmetric'),psfEff,'valid');
%PSF is always odd-sized, so we have to trim again
simIm  = simIm(2:end-1,2:end-1,2:end-1);

if showPlots
    convFig = figure;    
    subplot(1,2,1);
    imagesc(squeeze(max(simIm,[],1))),axis image,colormap gray
    title('Convolved Image, yz')
    subplot(1,2,2);
    imagesc(squeeze(max(simIm,[],3))),axis image,colormap gray    
    title('Convolved Image, xy')
end



 
   % ------------- Noise and Backgrond Generation --------- %


%Create the background noise and offset
imNoise = uint16(poissrnd(noiseLambda,p.ImageSize) + backOffset);

simIm = simIm + imNoise;

if showPlots
    noiseFig = figure;    
    subplot(2,2,1);
    imagesc(squeeze(max(simIm,[],1))),axis image,colormap gray
    title('Noised Image, yz MIP')
    subplot(2,2,2);
    imagesc(squeeze(max(simIm,[],3))),axis image,colormap gray    
    title('Noised Image, xy MIP')
    subplot(2,2,3);
    imagesc(squeeze(simIm(:,:,p.ImageSize(3)/2))),axis image,colormap gray    
    title('Noised Image, xy Slice')
    subplot(2,2,4);
    imagesc(squeeze(simIm(p.ImageSize(1)/2,:,:))),axis image,colormap gray    
    title('Noised Image, yz Slice')
end





