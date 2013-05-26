function [meanDispError,meanForceError,bead_x, bead_y, Av] = testSingleForce(f,d,minCorLength,dataPath,bead_x, bead_y, Av)
%% single force experiment
% input parameters to be replaced with function inputs
% f=2000; %Pa
% d=10;
% minCorLength = 21;
imgPath=[dataPath filesep 'Beads'];
refPath=[dataPath filesep 'Reference'];
orgPath=[dataPath filesep 'Original'];
analysisFolder = dataPath;
if ~exist(refPath,'dir') || ~exist(orgPath,'dir')
    mkdir(imgPath);
    mkdir(refPath);
    mkdir(orgPath);
end
%% reference image (300x200)
xmax=200;
ymax=300;
nPoints = 2000;
bead_r = 100; % nm
pixSize = 108; % nm/pix 60x
sigma = 2*bead_r/pixSize;
if nargin ==7
    refimg = simGaussianBeads(xmax,ymax, sigma, ...
        'x',bead_x,'y',bead_y,'A',Av, 'Border', 'truncated');
elseif nargin<5
    [refimg,bead_x, bead_y, ~, Av] = simGaussianBeads(xmax,ymax, sigma, ...
        'npoints', nPoints, 'Border', 'truncated');
end

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
% figure, imshow(refimg,[])
% figure, imshow(beadimg,[])

%% original force for comparison
force_x = assumedForceAniso2D(1,x_mat_u,y_mat_u,(100),150,0,f,d,d,forceType);
force_y = assumedForceAniso2D(2,x_mat_u,y_mat_u,(100),150,0,f,d,d,forceType);
%% saving original displacement field and force field
save([orgPath filesep 'data.mat'],'ux','uy','x_mat_u','y_mat_u','bead_x','bead_ux','bead_y','bead_uy','force_x','force_y');
%% save images
imwrite(uint16(refimg*2^16/max(max(refimg))),[refPath filesep 'img1ref.tif'],'tif')
imwrite(uint16(beadimg*2^16/max(max(beadimg))),[imgPath filesep 'img2bead.tif'],'tif')

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
%% Run the displacement field tracking
MD.getPackage(iPack).getProcess(2).run();
%% Create third process and run
MD.getPackage(iPack).createDefaultProcess(3)
params = MD.getPackage(iPack).getProcess(3).funParams_;
MD.getPackage(iPack).getProcess(3).setPara(params);
MD.getPackage(iPack).getProcess(3).run();

%% Create force reconstruction process and run
MD.getPackage(iPack).createDefaultProcess(4)
params = MD.getPackage(iPack).getProcess(4).funParams_;

params.YoungModulus = 8000;
params.regParam = 1e-6;
params.solMethodBEM = 'QR';
params.basisClassTblPath = '/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/TFM Basis Function SFT.mat';
MD.getPackage(iPack).getProcess(4).setPara(params);
MD.getPackage(iPack).getProcess(4).run();

MD.save;

%% Postprocessing - saving and analyzing force field
% Loading displacement field and force field
% Load the displField
displField=MD.getPackage(iPack).getProcess(3).loadChannelOutput;
% finding displacement at bead location
org_ux = zeros(size(displField(1).pos(:,1)));
org_uy = zeros(size(displField(1).pos(:,1)));
nmPoints = length(displField(1).pos(:,1));
for k=1:nmPoints
    [~,indcol_closest_x] = min(abs(x_mat_u(1,:)-displField(1).pos(k,1)),[],2);
    [~,indrow_closest_y] = min(abs(y_mat_u(:,1)-displField(1).pos(k,2)),[],1);
    row_bottom = max(1,indrow_closest_y-3);
    row_top = min(size(x_mat_u,1),indrow_closest_y+3);
    col_bottom = max(1,indcol_closest_x-3);
    col_top = min(size(y_mat_u,2),indcol_closest_x+3);
    loc_xmat = x_mat_u(row_bottom:row_top,col_bottom:col_top);
    loc_ymat = y_mat_u(row_bottom:row_top,col_bottom:col_top);
    loc_ux = ux(row_bottom:row_top,col_bottom:col_top);
    loc_uy = uy(row_bottom:row_top,col_bottom:col_top);
    org_ux(k) = interp2(loc_xmat,loc_ymat,loc_ux,displField(1).pos(k,1),displField(1).pos(k,2));
    if isnan(org_ux(k))
        org_ux(k) = ux(indrow_closest_y,indcol_closest_x);
    end
    org_uy(k) = interp2(loc_xmat,loc_ymat,loc_uy,displField(1).pos(k,1),displField(1).pos(k,2));
    if isnan(org_uy(k))
        org_uy(k) = uy(indrow_closest_y,indcol_closest_x);
    end
end
% errors in displacementfield
displField.vec(isnan(displField.vec(:,1)),:) = 0;
meanDispError= nansum(((org_ux-displField(1).vec(:,1)).^2+(org_uy-displField(1).vec(:,2)).^2).^.5); %is there a way to normalize this? by the number of beads?

% Load the forcefield
forceField=MD.getPackage(iPack).getProcess(4).loadChannelOutput;

% finding force at mesh location
org_fx = zeros(size(forceField(1).pos(:,1)));
org_fy = zeros(size(forceField(1).pos(:,1)));
nmfPoints = length(forceField(1).pos(:,1));
for k=1:nmfPoints
    [~,indcol_closest_x] = min(abs(x_mat_u(1,:)-forceField(1).pos(k,1)),[],2);
    [~,indrow_closest_y] = min(abs(y_mat_u(:,1)-forceField(1).pos(k,2)),[],1);
    row_bottom = max(1,indrow_closest_y-3);
    row_top = min(size(x_mat_u,1),indrow_closest_y+3);
    col_bottom = max(1,indcol_closest_x-3);
    col_top = min(size(y_mat_u,2),indcol_closest_x+3);
    loc_xmat = x_mat_u(row_bottom:row_top,col_bottom:col_top);
    loc_ymat = y_mat_u(row_bottom:row_top,col_bottom:col_top);
    loc_fx = force_x(row_bottom:row_top,col_bottom:col_top);
    loc_fy = force_y(row_bottom:row_top,col_bottom:col_top);
    org_fx(k) = interp2(loc_xmat,loc_ymat,loc_fx,forceField(1).pos(k,1),forceField(1).pos(k,2));
    if isnan(org_fx(k))
        org_fx(k) = force_x(indrow_closest_y,indcol_closest_x);
    end
    org_fy(k) = interp2(loc_xmat,loc_ymat,loc_fy,forceField(1).pos(k,1),forceField(1).pos(k,2));
    if isnan(org_fy(k))
        org_fy(k) = force_y(indrow_closest_y,indcol_closest_x);
    end
end

% errors in force field
meanForceError=nansum(((org_fx-forceField(1).vec(:,1)).^2+(org_fy-forceField(1).vec(:,2)).^2).^.5);
% heatmap creation and saving - i'll do it later

% 

return

%% input parameters to be replaced with function inputs
% f=2000; %Pa
% d=10;
% minCorLength = 21;
% dataPath='/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting';
% 
% testSingleForce(f,d,minCorLength,dataPath)
