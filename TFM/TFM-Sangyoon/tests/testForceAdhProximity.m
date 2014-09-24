function [detected,nDetected,fm1,fm2,fMap,cropInfo,bead_x, bead_y, Av] = testForceAdhProximity(d,f1,f2,r1,r2,method,dataPath, bead_x, bead_y, Av)
% testForceAdhProximity is a function that tests how  forces from two close
% adhesions are identified independently.
% input: 
%               d:              distance between adhesions (preferably even
%                                number)
%               f1:             force mag in adhesion 1 at left
%               f2:             force mag in adhesion 2 at right
%               r1:             adhesion radius in adhesion 1 at left
%               r2:             adhesion radius in adhesion 2 at right
%               method:     'L1' or 'L2'
%               dataPath:   data path to store all the results
% output: 
%               detected:  true if the two adhesion forces are identified
%                                 properly (adhesions close enough to force local maxima)
%               fm1:          force mag in adhesion 1 within the mesh element
%               fm2:          force mag in adhesion 2 within the mesh element
%               fMap:        traction Image

%% Preparing synthetic bead images
% reference image (200x200)
xmax=100;
ymax=100;
nPoints = 3000; % was 7000
% bead_r = 40; % nm
pixSize = 72; % nm/pix 90x
sigma = 1.68; % after getGaussianPSFsigma(NA,M,pixSize,lambda); 
if nargin<8
    bead_x=[];
end
if ~isempty(bead_x)
        refimg = simGaussianBeads(xmax,ymax, sigma, ...
        'x',bead_x,'y',bead_y,'A',Av, 'Border', 'truncated');
else
    Aorg = [300+100*randn(1,nPoints*4/5) 300+600*randn(1,nPoints*1/5)];
    Aorg(Aorg<0)=-Aorg(Aorg<0)+50;

    [refimg,bead_x, bead_y, ~, Av] = simGaussianBeads(xmax,ymax, sigma, ...
        'npoints', nPoints, 'Border', 'truncated','A',Aorg);
end

%% Noise addition (10%) % it was 5%
noiseLevel = 0.05;
refimg2 = 700+100*noiseLevel*randn(ymax,xmax) + refimg;% + 0.05*imgRange*(0.5-rand(ymax,xmax))*max(refimg(:));
% figure, imshow(refimg2,[])

% bead images
%% displacement field
E=8000;  %Young's modulus, unit: Pa
meshPtsFwdSol=2^10;
forceType = 'groupForce';

gridSpacing = 1;
xmin = gridSpacing;
ymin = gridSpacing;

[x_mat_u, y_mat_u]=meshgrid(xmin:gridSpacing:xmax,ymin:gridSpacing:ymax);

%% displacement field
midx=ceil(xmax/2);
midy=ceil(ymax/2);
posNA = [midx-d/2 midy;
                  midx+d/2 midy];

[ux, uy]=fwdSolution(x_mat_u,y_mat_u,E,xmin,xmax,ymin,ymax,...
    @(x,y) assumedForceAniso2D(1,x,y,midx-d/2,midy,f1,0,r1,r1,forceType)+...
    assumedForceAniso2D(1,x,y,midx+d/2,midy,f2,0,r2,r2,forceType),...
    @(x,y) assumedForceAniso2D(2,x,y,midx-d/2,midy,f1,0,r1,r1,forceType)+...
    assumedForceAniso2D(2,x,y,midx+d/2,midy,f2,0,r2,r2,forceType),'fft',[],meshPtsFwdSol);

force_x = assumedForceAniso2D(1,x_mat_u,y_mat_u,midx-d/2,midy,f1,0,r1,r1,forceType)+...
    assumedForceAniso2D(1,x_mat_u,y_mat_u,midx+d/2,midy,f2,0,r2,r2,forceType);
force_y = assumedForceAniso2D(2,x_mat_u,y_mat_u,midx-d/2,midy,f1,0,r1,r1,forceType)+...
    assumedForceAniso2D(2,x_mat_u,y_mat_u,midx+d/2,midy,f2,0,r2,r2,forceType);
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
beadimg = 700+100*noiseLevel*randn(ymax,xmax) + beadimg;% + 0.05*imgRange*(0.5-rand(ymax,xmax))*max(refimg(:));
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
save([orgPath filesep 'data.mat']);
%% save images
imwrite(uint16(refimg2*2^16/max(max(refimg2))),[refPath filesep 'img1ref.tif'],'tif')
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
params.highRes = true;
params.useGrid = false;
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

params.YoungModulus = 8000;
if strcmp(method,'L1')
    params.solMethodBEM = '1NormReg';
    params.regParam = 1.8e-5;
elseif strcmp(method,'L2')
    params.solMethodBEM = 'QR';
    params.regParam = 2.4e-7;
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
% forceMap
% force peak ratio
% maskForce = ((x_mat_u-100).^2+(y_mat_u-100).^2).^0.5<=d/2;
% forceForceIdx = maskVectors(forceField(1).pos(:,1),forceField(1).pos(:,2),maskForce);
% make  a an interpolated TF image and get the peak force because force
% mesh is sparse
[fMap,XI,YI]=generateHeatmapFromField(forceField);
% local maxima quantification
x1=posNA(1,1);
x2=posNA(2,1);
y1=posNA(1,2);
y2=posNA(2,2);
distThres = forceField.pos(2,2)-forceField.pos(2,1);
ynmin = round(y1)-YI(1,1)-distThres+1;
ynmax = round(y1)-YI(1,1)+distThres+1;
xnmin = round(x1)-XI(1,1)-distThres+1;
xnmax = round(x1)-XI(1,1)+distThres+1;
forceNeigh = fMap(ynmin:ynmax,xnmin:xnmax);
fm1 = max(forceNeigh(:));    
ynmin = round(y2)-YI(1,1)-distThres+1;
ynmax = round(y2)-YI(1,1)+distThres+1;
xnmin = round(x2)-XI(1,1)-distThres+1;
xnmax = round(x2)-XI(1,1)+distThres+1;
forceNeigh = fMap(ynmin:ynmax,xnmin:xnmax);
fm2 = max(forceNeigh(:));    
cropInfo = [XI(1,1), YI(1,1)];

% %new mask with XI and YI
% maskForceXIYI = ((XI-x1).^2+(YI-y1).^2).^0.5<=20 | ((XI-x2).^2+(YI-y2).^2).^0.5<=20;
% fImg = locmax2d(fMap,[d+20 20]);
%     % Identify adhesion location ([x1, y1], [x2, y2])
%     % See if they are close to force loc max
%     % if sum(sum((fImg>0).*maskForceXIYI))==2
% fImgCrop=(fImg>0).*maskForceXIYI;
% locmaxIdx = find(fImgCrop);
% if ~isempty(locmaxIdx)
%     xm = XI(locmaxIdx);
%     ym = YI(locmaxIdx);
%     [~,dist1] = KDTreeClosestPoint([xm ym], [x1 y1]);
%     [~,dist2] = KDTreeClosestPoint([xm ym], [x2 y2]);
% 
%     if dist1<7 && dist2<7
%         detected = 1;
%     elseif dist1<7 || dist2<7
%         detected = 0.5;
%     else
%         detected = 0;
%     end
% else
%     detected = 0;
% end
maskForceXIYI = ((XI-x1).^2+(YI-y1).^2).^0.5<=distThres | ((XI-x2).^2+(YI-y2).^2).^0.5<=distThres;
pstruct = pointSourceDetection(fMap,1.1,'mask',maskForceXIYI);
if ~isempty(pstruct)
    fMaxima = [pstruct.x'+cropInfo(1)-1 pstruct.y'+cropInfo(2)-1];
    nDetected = length(pstruct.x);
    [~,dist1] = KDTreeClosestPoint(fMaxima, [midx-d/2 midy]);
    [~,dist2] = KDTreeClosestPoint(fMaxima, [midx+d/2 midy]);
    if dist1<distThres && dist2<distThres
        detected = 1;
    elseif dist1<distThres || dist2<distThres
        detected = 0.5;
    else
        detected = 0;
    end
else
    detected = 0;
    nDetected = 0;
end

%% display original forcefield
h2 = generateHeatmapFromGridData(x_mat_u,y_mat_u,force_x,force_y,[dataPath '/Original forcefield'],0,400,false,size(XI,2),size(YI,1));
close(h2)

% if ((x1-x1m)^2+(y1-y1m)^2)^0.5<7 && ((x2-x2m)^2+(y2-y2m)^2)^0.5<7
%     detected = true;
% else
%     detected = false;
% end
% % else
% detected = false;
% end
return
