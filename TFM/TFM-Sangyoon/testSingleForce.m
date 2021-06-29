function [meanDispErrorAdh,meanDispErrorBG,dispDetec,meanForceErrorAdh,...
    meanForceErrorBG,peakForceRatio,forceDetec,beadsOnAdh,...
    bead_x, bead_y, Av] = testSingleForce(f,d,minCorLength,dataPath,...
                                        bead_x, bead_y, Av,varargin)
% Input check
ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('solMethodBEM','QR', @ischar);
ip.addParameter('regParam',1e-6,@isnumeric);
ip.addParameter('addNoise',false,@islogical);
ip.addParameter('whiteNoise',0.05,@isscalar);
ip.addParameter('nPoints',7000,@isscalar);
ip.addParameter('pathBasisClassTbl',[],@(x) isempty(x) || ischar(x));
ip.addParameter('E',8000,@isscalar);

ip.parse(varargin{:});
solMethodBEM = ip.Results.solMethodBEM;    
regParam = ip.Results.regParam;    
addNoise = ip.Results.addNoise;    
whiteNoise = ip.Results.whiteNoise;    
nPoints = ip.Results.nPoints;
pathBasisClassTbl = ip.Results.pathBasisClassTbl;
E = ip.Results.E; % in Pa

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
meshPtsFwdSol=2^8;
xmax=meshPtsFwdSol;
ymax=meshPtsFwdSol;

bead_r = 40; % nm
pixSize = 6.5e-6; %m/pix     %72e-9; % nm/pix 90x
NA = 1.49; %TIRF
lambda = 665e-9;
M = 60; %60 %1 because I use alreay magnified pixSize.
sigma = getGaussianPSFsigma(NA,M,pixSize,lambda); %1.68; % after 

if ~isempty(bead_x)
%     refimg = simGaussianSpots(xmax,ymax,sigma, ...
%         'x',bead_x,'y',bead_y,'A',Av, 'Border', 'truncated');
    refimg = simGaussianBeads(xmax,ymax, sigma, ...
        'x',bead_x,'y',bead_y,'A',Av, 'Border', 'truncated');
else
%     [refimg,bead_x, bead_y, ~, Av] = simGaussianSpots(xmax,ymax, sigma, ...
%         'npoints', nPoints, 'Border', 'truncated','A',0.3+rand(1,nPoints));
    curA = zeros(1,nPoints);
    kon=100; koff=50;
    nDye=650; dyeAmp=1/nDye;
    for pp=1:nPoints
        curA(1,pp) = dyeAmp * sum(rand(1,nDye)<kon/(kon+koff));
    end
    [refimg,bead_x, bead_y, ~, Av] = simGaussianBeads(xmax,ymax, sigma, ...
        'npoints', nPoints, 'Border', 'truncated','A',curA);
end
nPoints = length(bead_x);

%% Noise addition (5%)
refimg = refimg+0.05*rand(ymax,xmax)*max(refimg(:));
% %% Noise addition (10%)
% refimg = refimg+0.10*rand(ymax,xmax)*max(refimg(:));
%% Now displacement field from given force
% E=8000;  %Young's modulus, unit: Pa
forceType = 'groupForce';

gridSpacing = 1;
xmin = gridSpacing;
ymin = gridSpacing;

[x_mat_u, y_mat_u]=meshgrid(xmin:gridSpacing:xmax,ymin:gridSpacing:ymax);

%% original force for comparison
force_x = assumedForceAniso2D(1,x_mat_u,y_mat_u,(xmax/2),ymax/2,0,f,d,d,forceType);
force_y = assumedForceAniso2D(2,x_mat_u,y_mat_u,(xmax/2),ymax/2,0,f,d,d,forceType);
%% force noise
if addNoise
    force_x = 100*rand(ymax,xmax)-50 + force_x;
    force_y = 100*rand(ymax,xmax)-50 + force_y;
    [ux, uy]=fwdSolution(x_mat_u,y_mat_u,E,xmin,xmax,ymin,ymax,...
        force_x,force_y,'fft',[],meshPtsFwdSol,50000,0.5,false,true); %,'conv',[],meshPtsFwdSol);
else
    [ux, uy]=fwdSolution(x_mat_u,y_mat_u,E,xmin,xmax,ymin,ymax,...
        @(x,y) assumedForceAniso2D(1,x,y,(xmax/2),ymax/2,0,f,d,d,forceType),...
        @(x,y) assumedForceAniso2D(2,x,y,(xmax/2),ymax/2,0,f,d,d,forceType),'fft',[],meshPtsFwdSol); %,'conv',[],meshPtsFwdSol);
end
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

% beadimg = simGaussianSpots(xmax,ymax, sigma, ...
%     'x',bead_x+bead_ux,'y',bead_y+bead_uy,'A',Av,'Border', 'truncated');
beadimg = simGaussianBeads(xmax,ymax, sigma, ...
    'x',bead_x+bead_ux,'y',bead_y+bead_uy,'A',Av,'Border', 'truncated');
% figure, imshow(refimg,[])
% figure, imshow(beadimg,[])
%% Noise addition (5%)
beadimg = beadimg+whiteNoise*rand(ymax,xmax)*max(beadimg(:));
% %% Noise addition (10%)
% beadimg = beadimg+0.10*rand(ymax,xmax)*max(beadimg(:));

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
MD.pixelSize_=pixSize/M*1e9; %nm/pix
MD.camBitdepth_=16;
MD.timeInterval_ = 5;
MD.notes_=['Created for single force test purposes with f=' num2str(f) ' and d=' num2str(d)]; 

% Run sanityCheck on MovieData. 
% Check image size and number of frames are consistent. 
% Save the movie if successfull
MD.sanityCheck;
% Save the movie
MD.save;

%% Load the movie
% clear MD
% MD=MovieData.load(fullfile(analysisFolder,'movieData.mat'));

%% Create TFM package and retrieve package index
MD.addPackage(TFMPackage(MD));
iPack=  MD.getPackageIndex('TFMPackage');

%% Create second process
MD.getPackage(iPack).createDefaultProcess(2)
params = MD.getPackage(iPack).getProcess(2).funParams_;

%% Parameters in displacement field tracking
refFullPath = [refPath filesep 'img1ref.tif'];

params.referenceFramePath = refFullPath;
if f<300
    params.maxFlowSpeed = 5;
elseif f<1000
    params.maxFlowSpeed = 10;
elseif f<2000
    params.maxFlowSpeed = 20;
else
    params.maxFlowSpeed = 40;
end    
params.alpha = 0.05;
params.minCorLength = minCorLength;
params.mode = 'accurate';
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

params.YoungModulus = E;
params.method = 'FTTC';
params.regParam = regParam;
params.solMethodBEM = solMethodBEM;
params.useLcurve = false;
% params.basisClassTblPath = '/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/TFM Basis Function SFT.mat';
% params.basisClassTblPath = '/hms/scratch1/sh268/TFM Basis Function/TFM Basis Function SFT.mat';
% params.basisClassTblPath = '/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM basis functions/TFM Basis Function SFT.mat';
params.basisClassTblPath = pathBasisClassTbl;
MD.getPackage(iPack).getProcess(4).setPara(params);
MD.getPackage(iPack).getProcess(4).run();

MD.save;

%% Postprocessing - saving and analyzing force field
% Loading displacement field and force field
% Load the displField
disp('Calculating displacement errors and force errors...')
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
%% errors in displacementfield
% displField.vec(isnan(displField.vec(:,1)),:) = 0;
maskForce2 = ((x_mat_u-xmax/2).^2+(y_mat_u-ymax/2).^2).^0.5<=d/2*6;
dispIdx = maskVectors(displField(1).pos(:,1),displField(1).pos(:,2),maskForce2);

meanDispErrorAdh= nansum(((org_ux(dispIdx)-displField(1).vec(dispIdx,1)).^2+(org_uy(dispIdx)-displField(1).vec(dispIdx,2)).^2).^.5)/sum(~isnan((displField(1).vec(dispIdx,2)))); %normalized by the number of beads
meanDispErrorBG= nansum(((org_ux(~dispIdx)-displField(1).vec(~dispIdx,1)).^2+(org_uy(~dispIdx)-displField(1).vec(~dispIdx,2)).^2).^.5)/sum(~isnan((displField(1).vec(~dispIdx,2)))); %normalized by the number of beads
% detectability (u at force application / u at background)
maskForce = ((x_mat_u-xmax/2).^2+(y_mat_u-ymax/2).^2).^0.5<=d/2*2;
dispDetecIdx = maskVectors(displField(1).pos(:,1),displField(1).pos(:,2),maskForce);
if isempty(dispDetecIdx)
    dispDetec = 0;
else
    displFieldForce = displField(1).vec(dispDetecIdx,:);
    displFieldMag = (displFieldForce(:,1).^2+displFieldForce(:,2).^2).^0.5;
    backgroundIdx = maskVectors(displField(1).pos(:,1),displField(1).pos(:,2),~bwmorph(maskForce,'dilate',floor(d/2)));
    displFieldBgd = displField(1).vec(backgroundIdx,:);
    displFieldBgdMag = (displFieldBgd(:,1).^2+displFieldBgd(:,2).^2).^0.5;
    if isempty(displFieldMag)
        dispDetec = 0;
    else
        % sort and take top 10% mags
        displFieldBgdMagsorted = sort(displFieldBgdMag);
        dispDetec = mean(displFieldMag)/mean(displFieldBgdMagsorted(1:floor(0.1*length(displFieldBgdMagsorted))));%1:10));%f
    end
end
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

% heatmap creation and saving - i'll do it later

% force peak ratio
maskForce = ((x_mat_u-xmax/2).^2+(y_mat_u-ymax/2).^2).^0.5<=d/2;
% forceForceIdx = maskVectors(forceField(1).pos(:,1),forceField(1).pos(:,2),maskForce);
% make  a an interpolated TF image and get the peak force because force
% mesh is sparse
[fMap,XI,YI]=generateHeatmapFromField(forceField);
%new mask with XI and YI
maskForceXIYI = ((XI-xmax/2).^2+(YI-ymax/2).^2).^0.5<=d/2;

% if isempty(forceForceIdx)
%     peakForceRatio = 0;
% else
x_vec = reshape(x_mat_u,[],1);
y_vec = reshape(y_mat_u,[],1);
force_x_vec = reshape(force_x,[],1);
force_y_vec = reshape(force_y,[],1);
%     forceFieldForce = forceField(1).vec(forceForceIdx,:);
%     forceFieldMag = (forceFieldForce(:,1).^2+forceFieldForce(:,2).^2).^0.5;
fMapFiltered = fMap.*maskForceXIYI;
forceFieldMag = fMapFiltered(fMapFiltered>0);
orgFieldForceIdx = maskVectors(x_vec,y_vec,maskForce);
orgFieldForceMag = (force_x_vec(orgFieldForceIdx).^2+force_y_vec(orgFieldForceIdx).^2).^0.5;

backgroundIdx = maskVectors(forceField(1).pos(:,1),forceField(1).pos(:,2),~bwmorph(maskForce,'dilate',10));
forceFieldBgd = forceField(1).vec(backgroundIdx,:);
forceFieldBgdMag = (forceFieldBgd(:,1).^2+forceFieldBgd(:,2).^2).^0.5;
if isempty(forceFieldMag)
    peakForceRatio = 0;
    forceDetec = 0;
else
    peakForceRatio = mean(forceFieldMag)/mean(orgFieldForceMag);
    forceFieldBgdMag = sort(forceFieldBgdMag,'descend');
%     forceDetec = mean(forceFieldMag)/mean(forceFieldBgdMag(1:round(length(forceFieldMag)/2)));
    forceDetec = max(forceFieldMag)/max(forceFieldBgdMag(1:length(forceFieldMag)));
end
%% errors in force field
forceIdx = maskVectors(forceField(1).pos(:,1),forceField(1).pos(:,2),maskForce);
meanForceErrorAdh=nansum(((org_fx(forceIdx)-forceField(1).vec(forceIdx,1)).^2+(org_fy(forceIdx)-forceField(1).vec(forceIdx,2)).^2).^.5)/sum(~isnan((forceField(1).vec(forceIdx,2))));
meanForceErrorBG=nansum(((org_fx(backgroundIdx)-forceField(1).vec(backgroundIdx,1)).^2+(org_fy(backgroundIdx)-forceField(1).vec(backgroundIdx,2)).^2).^.5)/sum(~isnan((forceField(1).vec(backgroundIdx,2))));
%% beadsOnAdh
beadIdx = maskVectors(displField(1).pos(:,1),displField(1).pos(:,2),maskForce);
if sum(beadIdx)
    beadsOnAdh = true;
else
    beadsOnAdh = false;
end
return

%% input parameters to be replaced with function inputs
% f=2000; %Pa
% d=10;
% minCorLength = 21;
% dataPath='/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting';
% 
% testSingleForce(f,d,minCorLength,dataPath)
