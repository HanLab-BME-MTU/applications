%% single force experiment
clear
close all

%% input parameters to be replaced with function inputs
f=2000; %Pa
d=10;
minCorLength = 21;
imgPath='/Users/joshua2/Documents/PostdocResearch/Traction Force/testSingleForce/f2000d10/Beads';
refPath='/Users/joshua2/Documents/PostdocResearch/Traction Force/testSingleForce/f2000d10/Reference';
analysisFolder = '/Users/joshua2/Documents/PostdocResearch/Traction Force/testSingleForce/f2000d10';
if ~exist(refPath,'dir')
    mkdir(imgPath);
    mkdir(refPath);
end
%% reference image (300x200)
xmax=200;
ymax=300;
nPoints = 2000;
bead_r = 60; % nm
pixSize = 108; % nm/pix 60x
sigma = 2*bead_r/pixSize;
[refimg,bead_x, bead_y, sv, Av] = simGaussianBeads(xmax,ymax, sigma, ...
    'npoints', nPoints, 'Border', 'truncated');

%% Now displacement field from given force
E=8000;  %Young's modulus, unit: Pa
meshPtsFwdSol=2^10;
forceType = 'groupForce';

gridSpacing = 1;
xmin = gridSpacing;
ymin = gridSpacing;

[x_mat_u, y_mat_u]=meshgrid(xmin:gridSpacing:xmax,ymin:gridSpacing:ymax);


[ux, uy]=fwdSolution(x_mat_u,y_mat_u,E,xmin,xmax,ymin,ymax,...
    @(x,y) assumedForceAniso2D(1,x,y,(100),150,0,f,d,d,forceType),...
    @(x,y) assumedForceAniso2D(2,x,y,(100),150,0,f,d,d,forceType),'fft',[],meshPtsFwdSol); %,'conv',[],meshPtsFwdSol);

% figure, quiver(x_mat_u,y_mat_u,ux,uy)

%% bead image
% finding displacement at bead location
bead_ux = zeros(size(bead_x));
bead_uy = zeros(size(bead_y));
for k=1:nPoints
    [~,indcol_closest_x] = min(abs(x_mat_u(1,:)-bead_x(k)),[],2);
    [~,indrow_closest_y] = min(abs(y_mat_u(:,1)-bead_y(k)),[],1);
    row_bottom = max(1,indrow_closest_y-2);
    row_top = min(size(x_mat_u,1),indrow_closest_y+2);
    col_bottom = max(1,indcol_closest_x-2);
    col_top = min(size(y_mat_u,2),indcol_closest_x+2);
    loc_xmat = x_mat_u(row_bottom:row_top,col_bottom:col_top);
    loc_ymat = y_mat_u(row_bottom:row_top,col_bottom:col_top);
    loc_ux = ux(row_bottom:row_top,col_bottom:col_top);
    loc_uy = uy(row_bottom:row_top,col_bottom:col_top);
    bead_ux(k) = interp2(loc_xmat,loc_ymat,loc_ux,bead_x(k),bead_y(k));
    if isnan(bead_ux(k))
        bead_ux(k) = ux(indrow_closest_y,indcol_closest_x);
    end
    bead_uy(k) = interp2(loc_xmat,loc_ymat,loc_uy,bead_x(k),bead_y(k));
    if isnan(bead_uy(k))
        bead_uy(k) = uy(indrow_closest_y,indcol_closest_x);
    end
end

beadimg = simGaussianBeads(xmax,ymax, sigma, ...
    'x',bead_x+bead_ux,'y',bead_y+bead_uy,'A',Av,'Border', 'truncated');
figure, imshow(refimg,[])
figure, imshow(beadimg,[])

imwrite(uint16(refimg*2^16/max(max(refimg))),[refPath filesep 'img1ref.tif'],'tif')
imwrite(uint16(beadimg*2^16/max(max(beadimg))),[imgPath filesep 'img2bead.tif'],'tif')

%% original force for comparison
force_x = assumedForceAniso2D(1,x_mat_u,y_mat_u,(100),150,0,f,d,d,forceType);
force_y = assumedForceAniso2D(2,x_mat_u,y_mat_u,(100),150,0,f,d,d,forceType);

%% Now force reconstruction via movieData (non-GUI mode)
% Retrieve current location
fullPath = [imgPath filesep 'img2bead.tif'];
path = fileparts(fullPath);
dataFolder=path;
% dataFolder=[imgPath filesep 'img2bead.tif']; %imgPath;
%% Channel creation
% Create a channels object
channel = Channel(dataFolder);
channel.fluorophore_='alexa647';
channel.emissionWavelength_=name2wavelength('alexa647')*1e9;
channel.imageType_='TIRF';
%% MovieData creation
% Constructor needs an array of channels and an output directory (for analysis)
MD = MovieData(channel,analysisFolder);

% Set the path where to store the MovieData object.
MD.setPath(analysisFolder);
MD.setFilename('movieData.mat');

% Set some additional movie properties
MD.numAperture_=1.49;
MD.pixelSize_=108;
MD.camBitdepth_=16;
MD.timeInterval_ = 5;
MD.notes_='Created for single force test purposes'; 

% Run sanityCheck on MovieData. 
% Check image size and number of frames are consistent. 
% Save the movie if successfull
MD.sanityCheck;
% Save the movie
MD.save;

%% Load the movie
clear MD
MD=MovieData.load(fullfile(analysisFolder,'movieData.mat'));

%% Create TFM package and retrieve package index
MD.addPackage(TFMPackage(MD));
iPack=  MD.getPackageIndex('TFMPackage');

%% Create second process
MD.getPackage(iPack).createDefaultProcess(2)
params = MD.getPackage(iPack).getProcess(2).funParams_;

%% Parameters in displacement field tracking
refFullPath = [refPath filesep 'img1ref.tif'];

params.referenceFramePath = refFullPath;
params.maxFlowSpeed = 40;
params.alpha = 0.05;
params.minCorLength = minCorLength;
MD.getPackage(iPack).getProcess(2).setPara(params);

MD.getPackage(iPack).getProcess(2).run();



