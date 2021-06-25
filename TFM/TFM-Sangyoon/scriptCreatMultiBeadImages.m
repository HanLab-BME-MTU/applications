%scriptCreateMultiBeadImages creates the same bead images that has
%fluctuating intensities. From Stubb et al:
dataPath = '/Volumes/GoogleDrive/My Drive/Documents/Grant/2021 NSF CAREER/Sci/Figs/beadDensitySimulation/imgStack';
mkdir(dataPath)
bead_d = 40; % nm
pixSize = 108e-9; % m/pix 90x
NA = 1.49; %TIRF
lambda = 665e-9;
M = 1; %1 because I use alreay magnified pixSize.

sigma = getGaussianPSFsigma(NA,M,pixSize,lambda);
% sigma = 1.146; % after getGaussianPSFsigma(NA,M,pixSize,lambda);

%% image size
% The positions of the beads were randomly distributed over the chosen
% field-of- view size (here 50 × 50 μm2). 
x_um = 30; y_um=30;
x_pix = round(x_um/pixSize*1e-6);
y_pix = round(y_um/pixSize*1e-6);
%% bead density
% The number of beads was chosen to obtain the desired density. 
beadDensityDesired = 5; %beads/um2
areaImg = x_um * y_um; % um2
nPoints = beadDensityDesired * areaImg;
%% image
% Each bead was simulated as a group of ∼650
% dyes homogeneously distributed in a 40 nm sphere. 
nDye = 650;
dyeAmp = 1/nDye;
% Each dye was allowed to blink independently with on/off rate of 100 s−1
% and 50 s−1 respectively over the entire acquisition without bleaching
% (200 frames at 10 ms exposure). 
kon = 100; %s-1
koff = 50; %s-1
curA = zeros(1,nPoints);
%% iteration
nImgs = 100;
padZeros=floor(log10(nImgs))+1;
for k=1:nImgs
    for ii=1:nPoints
        curA(1,ii) = dyeAmp * sum(rand(1,nDye)<kon/(kon+koff));
    end
    % The simulated fluorescence image produced
    % by this distribution of beads was created by convolution with a Gaussian
    % kernel with σ = 0.21 × λ ,54 where λ is the wavelength of emission NA
    % (here 700 nm) and NA is the numerical aperture of the microscope 
    if k==1
        [img,bead_x, bead_y, ~, ~] = simGaussianBeads(x_pix, y_pix, sigma, ...
            'npoints', nPoints, 'Border', 'truncated','A',curA);
    else
        img = simGaussianBeads(x_pix, y_pix, sigma, ...
            'x',bead_x,'y',bead_y,'A', curA, 'Border', 'truncated');
    end
    % save img
    imwrite(uint16(img*2^12), [dataPath filesep 'img' num2str(k,['%0.',int2str(padZeros),'d']) '.tif'],'tiff','Compression','none')
end
