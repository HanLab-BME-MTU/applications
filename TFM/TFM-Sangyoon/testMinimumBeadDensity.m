% testMinimumBeadDensity calculates force detectability depending on bead
% density and white noise.
clear
%% Simulations - with pixel-wise tracking
numBeads = [250 500 750 1000 1500 2000 2500 3000 4000];
nNB = length(numBeads);
pixSize = 108; %nm/px
pixSizeUm = pixSize/1000;
d = 4; %pix (which means ~400 nm)
noiseLevels = 0:0.05:0.2;
nNL = length(noiseLevels);
nExp = 10;
E=2000; %Pa

d_err_Adh = zeros(nNB,nNL,nExp);
d_err_BG = zeros(nNB,nNL,nExp);
f_err_ADh = zeros(nNB,nNL,nExp);
f_err_BG = zeros(nNB,nNL,nExp);
dispDetec = zeros(nNB,nNL,nExp);
forceDetec = zeros(nNB,nNL,nExp);
pFR = zeros(nNB,nNL,nExp);
beadsOnAdh = zeros(nNB, nExp);
cL = 15;%[9 15 21]
adhArea = pi/4*(d*pixSize)^2*1e-18; %m2
f = (15e-12)/adhArea; %15 picoNewton/adhesion area, Pa
adhAreaUm2 = adhArea*1e12;
w = 2^8; h=w;
fieldArea = w*h; %pix
fieldAreaUm2 = fieldArea*pixSizeUm^2;
beadDensity = numBeads/fieldAreaUm2; % #/um2
B2Bspacing = 1./beadDensity.^2;
% kk=0;
%% simulation for numBeads and noiseLevel
% pathBasisClassTbl='/storage/network/TFM_Development/TFM2D/TFM_basisClass/basisClass2kPa9pix.mat';
% pathBasisClassTbl='/storage/network/TFM_Development/TFM2D/TFM_basisClass/basisClass2kPa9pix.mat';
pathBasisClassTbl='$PROJECT/TFM/TFM_basisClass/basisClass2kPa9pix.mat';

for epm=1:nExp
    p=0;
    ii=0;
    for nB=numBeads
        ii=ii+1;
        jj=0;
        for curWN = noiseLevels
            jj=jj+1;
            p=p+1;
%             dataPath=['/storage/network/TFM_Development/TFM2D/BeadDensityRequirement/2kPaSimulation/f' num2str(round(f)) ...
            dataPath=['$PROJECT/TFM/TFMBeadDensity/2kPaSimulation/f' num2str(round(f)) ...
                'd' num2str(d) 'nB' num2str(nB) 'noise' num2str(curWN) 'exp' num2str(epm)];

            if jj==1
                [d_err_Adh(ii,jj,epm),d_err_BG(ii,jj,epm),dispDetec(ii,jj,epm),...
                    f_err_ADh(ii,jj,epm),f_err_BG(ii,jj,epm),pFR(ii,jj,epm),forceDetec(ii,jj,epm),...
                    beadsOnAdh(ii,epm),bead_x{ii,epm}, bead_y{ii,epm}, Av{ii,epm}]= ...
                    testSingleForce(f,d,cL,dataPath,[],[],[],'QR','whiteNoise',curWN,'nPoints',nB,...
                    'pathBasisClassTbl',pathBasisClassTbl,'E',E,'regParam',4e-4);
            else
                [d_err_Adh(ii,jj,epm),d_err_BG(ii,jj,epm),dispDetec(ii,jj,epm),...
                    f_err_ADh(ii,jj,epm),f_err_BG(ii,jj,epm),pFR(ii,jj,epm),forceDetec(ii,jj,epm),...
                    beadsOnAdh(ii,epm)]= ...
                    testSingleForce(f,d,cL,dataPath,bead_x{ii,epm}, bead_y{ii,epm}, Av{ii,epm},...
                    'QR','whiteNoise',curWN, 'pathBasisClassTbl',pathBasisClassTbl,'E',E,'regParam',4e-4);
            end
        end
    end
end
%% plotting
myCs=distinguishable_colors(nNL);
figure, hold on
for jj=1:nNL
    curCol = myCs(jj,:);
    brighterCol = curCol.*1.5;
    brighterCol(brighterCol>1) = 1;
    plot(beadDensity,squeeze(forceDetec(:,jj,:)),'Color',brighterCol);
    h(jj) = plot(beadDensity,mean(forceDetec(:,jj,:),3),'LineWidth',2,'Color',myCs(jj,:));
end
legend(h,'0 %','5 %','10%','15%','20%')
savefig('$PROJECT/TFM/TFMBeadDensity/2kPaSimulation/beadDensityVsForceDetec.fig')
% savefig('/storage/network/TFM_Development/TFM2D/BeadDensityRequirement/2kPaSimulation/beadDensityVsForceDetec.fig')
%% save
% dataPath='/storage/network/TFM_Development/TFM2D/BeadDensityRequirement/2kPaSimulation/allData.mat';
dataPath='$PROJECT/TFM/TFMBeadDensity/2kPaSimulation/allData.mat';
save(dataPath)
%% Sample multi-image creation
% I chose to use jj=2, epm 5, nB=4000 as representative example to apply
% LiveSRRF, for which I need to generate 100 images with some fluctuation.
% First, ref images:
if 0
    nImgs=100;
    jj=3; epm=5; ii=9;
    nB=numBeads(ii);
    kon=100; koff=50;
    nDye=650; dyeAmp=1/nDye;
    curWN=noiseLevels(jj);
    meshPtsFwdSol=2^8;
    w=meshPtsFwdSol;
    h=meshPtsFwdSol;
    bead_d = 40; % nm
    pixSize = 108e-9; % m/pix 60x
    NA = 1.49; %TIRF
    lambda = 665e-9;
    M = 1; %1 because I use alreay magnified pixSize.
    pixSizeNm=pixSize*1e9;

    sigma = getGaussianPSFsigma(NA,M,pixSize,lambda);

    dataPath=['/storage/network/TFM_Development/TFM2D/BeadDensityRequirement/2kPaSimulation/f' ...
        num2str(round(f)) 'd' num2str(d) 'nB' num2str(nB) 'noise' num2str(curWN) ...
        'exp' num2str(epm) filesep 'MultiImages' filesep 'Ref'];
    mkdir(dataPath)
    padZeros=floor(log10(nImgs))+1;

    for k=1:nImgs
        curA=zeros(1,nB);
        for pp=1:nB
            curA(1,pp) = dyeAmp * sum(rand(1,nDye)<kon/(kon+koff));
        end
        % The simulated fluorescence image produced
        % by this distribution of beads was created by convolution with a Gaussian
        % kernel with ? = 0.21 × ? ,54 where ? is the wavelength of emission NA
        % (here 700 nm) and NA is the numerical aperture of the microscope 
        img = simGaussianBeads(w, h, sigma, ...
            'x',bead_x{ii,epm},'y',bead_y{ii,epm},'A', curA, 'Border', 'truncated');
        % camera noise addition
        img = img+curWN*rand(h,w)*max(img(:));

        % save img
        imwrite(uint16(img*2^12), [dataPath filesep 'img' num2str(k,['%0.',int2str(padZeros),'d']) '.tif'],'tiff','Compression','none')
    end
    %% Now deformed bead iamge
    % load bead_ux and bead_uy
    dataPath=['/storage/network/TFM_Development/TFM2D/BeadDensityRequirement/2kPaSimulation/f' ...
        num2str(round(f)) 'd' num2str(d) 'nB' num2str(nB) 'noise' num2str(curWN) ...
        'exp' num2str(epm) filesep 'Original'];
    orgData = load([dataPath filesep 'data.mat']);

    dataPath=['/storage/network/TFM_Development/TFM2D/BeadDensityRequirement/2kPaSimulation/f' ...
        num2str(round(f)) 'd' num2str(d) 'nB' num2str(nB) 'noise' num2str(curWN) ...
        'exp' num2str(epm) filesep 'MultiImages' filesep 'BeadMulti'];
    mkdir(dataPath)

    for k=1:nImgs
        curA=zeros(1,nB);
        for pp=1:nB
            curA(1,pp) = dyeAmp * sum(rand(1,nDye)<kon/(kon+koff));
        end
        % The simulated fluorescence image produced
        % by this distribution of beads was created by convolution with a Gaussian
        % kernel with ? = 0.21 × ? ,54 where ? is the wavelength of emission NA
        % (here 700 nm) and NA is the numerical aperture of the microscope 
        img = simGaussianBeads(w, h, sigma, ...
            'x',bead_x{ii,epm}+orgData.bead_ux,'y',bead_y{ii,epm}+orgData.bead_uy,...
            'A', curA, 'Border', 'truncated');
        % camera noise addition
        img = img+curWN*rand(h,w)*max(img(:));

        % save img
        imwrite(uint16(img*2^12), [dataPath filesep 'img' num2str(k,['%0.',int2str(padZeros),'d']) '.tif'],'tiff','Compression','none')
    end
    %% Channel creation
    % Create a channels object
    dataPath=['/storage/network/TFM_Development/TFM2D/BeadDensityRequirement/2kPaSimulation/f' ...
        num2str(round(f)) 'd' num2str(d) 'nB' num2str(nB) 'noise' num2str(curWN) ...
        'exp' num2str(epm) filesep 'MultiImages' filesep 'Bead'];

    channel = Channel(dataPath);
    channel.fluorophore_='alexa647';
    channel.emissionWavelength_=name2wavelength('alexa647')*1e9;
    channel.imageType_='TIRF';
    %% MovieData creation
    % Constructor needs an array of channels and an output directory (for analysis)
    analysisFolder=['/storage/network/TFM_Development/TFM2D/BeadDensityRequirement/2kPaSimulation/f' ...
        num2str(round(f)) 'd' num2str(d) 'nB' num2str(nB) 'noise' num2str(curWN) ...
        'exp' num2str(epm) filesep 'MultiImages'];
    MD = MovieData(channel,analysisFolder);

    % Set the path where to store the MovieData object.
    MD.setPath(analysisFolder);
    MD.setFilename('movieData.mat');

    % SR mag=5
    magSR=5;
    % Set some additional movie properties
    MD.numAperture_=1.49;
    MD.pixelSize_=108/magSR;
    MD.camBitdepth_=16;
    MD.timeInterval_ = 5;
    MD.notes_=['Created for single force test purposes with LiveSRRF']; 

    % Run sanityCheck on MovieData. 
    % Check image size and number of frames are consistent. 
    % Save the movie if successfull
    MD.sanityCheck;
    % Save the movie
    MD.save;
    %% Create TFM package and retrieve package index
    MD.addPackage(TFMPackage(MD));
    iPack=  MD.getPackageIndex('TFMPackage');

    %% Create second process
    MD.getPackage(iPack).createDefaultProcess(2)
    params = MD.getPackage(iPack).getProcess(2).funParams_;

    %% Parameters in displacement field tracking
    refPath=['/storage/network/TFM_Development/TFM2D/BeadDensityRequirement/2kPaSimulation/f' ...
        num2str(round(f)) 'd' num2str(d) 'nB' num2str(nB) 'noise' num2str(curWN) ...
        'exp' num2str(epm) filesep 'MultiImages' filesep 'Ref'];
    refFullPath = [refPath ' - liveSRRF (AVG).tif'];

    params.referenceFramePath = refFullPath;
    params.maxFlowSpeed = 1*magSR;
    minCorLength = 21;
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

    %% Create force reconstruction process 
    MD.getPackage(iPack).createDefaultProcess(4)
    %% and run
    params = MD.getPackage(iPack).getProcess(4).funParams_;
    regParam = 1e-4;
    solMethodBEM = 'backslash';
    params.method = 'FTTC';
    params.YoungModulus = E;
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
    x_mat_u=orgData.x_mat_u;
    y_mat_u=orgData.y_mat_u;
    xmax=w; ymax=h;
    disp('Calculating displacement errors and force errors...')
    displField=MD.getPackage(iPack).getProcess(3).loadChannelOutput;
    % Load the forcefield
    forceField=MD.getPackage(iPack).getProcess(4).loadChannelOutput;
    % force peak ratio
    maskForce = ((x_mat_u-xmax/2).^2+(y_mat_u-ymax/2).^2).^0.5<=d/2;
    % forceForceIdx = maskVectors(forceField(1).pos(:,1),forceField(1).pos(:,2),maskForce);
    % make  a an interpolated TF image and get the peak force because force
    % mesh is sparse
    [fMap,XI,YI]=generateHeatmapFromField(forceField);
    % need to downscale fMap with magSR

    %new mask with XI and YI
    maskForceXIYI = ((XI-xmax/2).^2+(YI-ymax/2).^2).^0.5<=d/2;

    % if isempty(forceForceIdx)
    %     peakForceRatio = 0;
    % else
    x_vec = reshape(x_mat_u,[],1);
    y_vec = reshape(y_mat_u,[],1);
    force_x_vec = reshape(orgData.force_x,[],1);
    force_y_vec = reshape(orgData.force_y,[],1);
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

end
