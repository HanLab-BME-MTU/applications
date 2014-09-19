function [detected,fm1,fm2,tImg] = testForceAdhProximity(d,f1,f2,r1,r2,method,dataPath)
% testForceAdhProximity is a function that tests how  forces from two close
% adhesions are identified independently.
% input: 
%               d:              distance between adhesions
%               f1:             force mag in adhesion 1 at left
%               f2:             force mag in adhesion 2 at right
%               r1:             adhesion radius in adhesion 1 at left
%               r2:             adhesion radius in adhesion 2 at right
%               method:   'L1' or 'L2'
%               dataPath:  data path to store all the results
% output: 
%               detected:  true if the two adhesion forces are identified
%                                 properly
%               fm1:          force mag in adhesion 1 within the mesh element
%               fm2:          force mag in adhesion 2 within the mesh element
%               tImg:         traction Image
%% Preparing synthetic bead images
% reference image (200x200)
xmax=200;
ymax=200;
nPoints = 4000; % was 25000
bead_r = 40; % nm
pixSize = 72; % nm/pix 90x
sigma = 1.68; % after getGaussianPSFsigma(NA,M,pixSize,lambda); 
Aorg = [300+100*randn(1,nPoints*4/5) 300+600*randn(1,nPoints*1/5)];
Aorg(Aorg<0)=-Aorg(Aorg<0)+50;
[refimg,bead_x, bead_y, ~, Av] = simGaussianBeads(xmax,ymax, sigma, ...
        'npoints', nPoints, 'Border', 'truncated','A',Aorg);
%% Noise addition (10%) % it was 5%
noiseLevel = 0.1;
refimg2 = 700+700*noiseLevel*randn(ymax,xmax) + refimg;% + 0.05*imgRange*(0.5-rand(ymax,xmax))*max(refimg(:));
figure, imshow(refimg2,[])

% bead images
%% displacement field
E=8000;  %Young's modulus, unit: Pa
meshPtsFwdSol=2^10;
forceType = 'groupForce';

gridSpacing = 1;
xmin = gridSpacing;
ymin = gridSpacing;

[x_mat_u, y_mat_u]=meshgrid(xmin:gridSpacing:xmax,ymin:gridSpacing:ymax);

%% temporary - get the coordinates
posNA = [101-d/2 101;
                  101+d/2 101];
%% displacement field

[ux, uy]=fwdSolution(x_mat_u,y_mat_u,E,xmin,xmax,ymin,ymax,...
    @(x,y) assumedForceAniso2D(1,x,y,101-d/2,101,f1,0,r1,r1,forceType)+...
    assumedForceAniso2D(1,x,y,101+d/2,101,f1,0,r2,r2,forceType),...
    @(x,y) assumedForceAniso2D(2,x,y,101-d/2,101,f1,0,r1,r1,forceType)+...
    assumedForceAniso2D(2,x,y,101+d/2,101,f1,0,r2,r2,forceType),'fft',[],meshPtsFwdSol);

%% finding displacement at bead location
nPoints = length(bead_x);
bead_ux = zeros(size(bead_x));
bead_uy = zeros(size(bead_y));
for k=1:nPoints
    [~,indcol_closest_x] = min(abs(x_mat_u(1,:)-bead_x(k)),[],2);
    [~,indrow_closest_y] = min(abs(y_mat_u(:,1)-bead_y(k)),[],1);
    row_bottom = max(1,indrow_closest_y-2);
    row_top = min(size(x_mat_u,2),indrow_closest_y+2);
    col_bottom = max(1,indcol_closest_x-2);
    col_top = min(size(y_mat_u,1),indcol_closest_x+2);
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

% pixelSize = 0.108; % assuming 60x objective um/pixel
beadimg = simGaussianBeads(xmax,ymax, sigma,'x',bead_x+bead_ux,'y',bead_y+bead_uy,'A',Av, 'Border', 'truncated');
%% Noise addition (10%) % it was 5% before
% beadimg = beadimg+0.1*rand(ymax,xmax)*max(beadimg(:));
beadimg = 700+700*noiseLevel*randn(ymax,xmax) + beadimg;% + 0.05*imgRange*(0.5-rand(ymax,xmax))*max(refimg(:));
%% saving
imgPath=[dataPath filesep 'Beads'];
refPath=[dataPath filesep 'Reference'];
orgPath=[dataPath filesep 'Original'];
analysisFolder = dataPath;
if ~exist(refPath,'dir') || ~exist(orgPath,'dir')
    mkdir(imgPath);
    mkdir(refPath);
    mkdir(orgPath);
end
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
MD.pixelSize_=pixSize;
MD.camBitdepth_=16;
MD.timeInterval_ = 5;
MD.notes_=['Created for single force test purposes with f1=' num2str(f1) ' and d=' num2str(d)]; 

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
params.maxFlowSpeed = 10;
params.alpha = 0.05;
params.minCorLength = 17;
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
if strcmp(method,'L1')
    params.solMethodBEM = '1NormReg';
    params.regParam = 1e-4;
elseif strcmp(method,'L2')
    params.solMethodBEM = 'QR';
    params.regParam = 1e-7;
else
    display('The method should be either L1 or L2. The input does not belong to any of those. We use L2 as a default.')
    params.solMethodBEM = 'QR';
end
params.method = 'FastBEM';
params.useLcurve = false;
params.basisClassTblPath = '/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM basis functions/basisClass8kPaSimul.mat';
MD.getPackage(iPack).getProcess(4).setPara(params);
MD.getPackage(iPack).getProcess(4).run();

MD.save;

%% Postprocessing - saving and analyzing force field
% Loading displacement field and force field
% Load the displField
disp('Detecting local maxima in reconstructed force ... ')
% Load the forcefield
forceField=MD.getPackage(iPack).getProcess(4).loadChannelOutput;

% force peak ratio
% maskForce = ((x_mat_u-100).^2+(y_mat_u-100).^2).^0.5<=d/2;
% forceForceIdx = maskVectors(forceField(1).pos(:,1),forceField(1).pos(:,2),maskForce);
% make  a an interpolated TF image and get the peak force because force
% mesh is sparse
[fMap,XI,YI]=generateHeatmapFromField(forceField);
%new mask with XI and YI
maskForceXIYI = ((XI-100).^2+(YI-150).^2).^0.5<=d/2;

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
    forceDetec = mean(forceFieldMag)/mean(forceFieldBgdMag(1:round(length(forceFieldMag)/2)));
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
