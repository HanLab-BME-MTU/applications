function [rmsError,EOA,lcurvePath,fMap,cropInfo,orgMap,bead_x, bead_y, Av] = testForceAdhDensity(n,method,dataPath,baseDataPath, bead_x, bead_y, Av)
% testForceAdhProximity is a function that tests how  forces from two close
% adhesions are identified independently.
% input: 
%               d:              distance between adhesions (preferably even
%                                number)
%               f1:             force mag in adhesion 1 at left
%               r1:             adhesion radius in adhesion 1 at left
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
xmax=200;
ymax=200;
nPoints = 12000; % was 7000
% bead_r = 40; % nm
pixSize = 72; % nm/pix 90x
sigma = 1.68; % after getGaussianPSFsigma(NA,M,pixSize,lambda); 
if nargin<4
    baseDataPath=[];
    bead_x=[];
elseif nargin<5
    bead_x=[];
end
%% Bead image creation
if isempty(baseDataPath) % in case you build images and reconstruct forces out of them.
    if ~isempty(bead_x)
            refimg = simGaussianBeads(xmax,ymax, sigma, ...
            'x',bead_x,'y',bead_y,'A',Av, 'Border', 'truncated');
    else
%         Aorg = [300+100*randn(1,nPoints*3/5) 300+600*rand(1,nPoints*2/5)];
%         Aorg(Aorg<300)=300+400*rand(1,sum(Aorg<300));
        Aorg = [300+600*rand(1,nPoints)];
%         Aorg(Aorg<300)=300+400*rand(1,sum(Aorg<300));

        [refimg,bead_x, bead_y, ~, Av] = simGaussianBeads(xmax,ymax, sigma, ...
            'npoints', nPoints, 'Border', 'truncated','A',Aorg);
    end

    %% Noise addition (10%) % it was 5%
    noiseLevel = 0.1;
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

%% position definition - square grid
%     posx_min = 4*rmax; posx_max = xmax-4*rmax;
%     posy_min = 4*rmax; posy_max = ymax-4*rmax;
%     posx_vec = posx_min:4*rmax:posx_max;
%     posy_vec = posy_min:4*rmax:posy_max;
%     [posx,posy] = meshgrid(posx_vec,posy_vec);
%% position definition - random
    % radius determination
    rmin = 4;
    rmax = 9;
    rx = [rmin+rmin*rand(1,round(n*2/3)), (rmax)*rand(1,round(n/3))];
    rx(rx<rmin) = rmin+2*rmin*rand(1,sum(rx<rmin));
    ry = rx/rmin.*rx +rmax*rand(1,n);
    rx(rx>0.9*rmax & ry./rx>4) = rx(rx>0.9*rmax & ry./rx>4)*0.8;
%     figure, hist(rx,20)
%     figure, hist(ry,20)
    % magnitude 
    fmin = 100; % Pa
    fmax= 1600;
%     f = [fmin+2*fmin*rand(1,round(n*1/3)), (fmax)*rand(1,round(n*2/3))];
    f = fmin+(fmax-fmin)*rand(1,n);
    f(f<fmin) = fmin +(fmin+fmax)/2*rand(1,sum(f<fmin));
    % Estimate coverage area with rx, ry.
    areaAdhesions = sum(pi*rx.*ry);
    if areaAdhesions>0.7*xmax*ymax
        display('n is too much for adhesion distribution. Please reduce n.')
    end
%% now placing the adhesion one by one
    distMin=10;
    posx = zeros(n,1);
    posy = zeros(n,1);
    iMax= 1000; % iteration maximum
    bandWith = 10;
    posDiscard = [];
    for ii=1:n
        posx_min = bandWith+rx(ii); posx_max = xmax-(bandWith+rx(ii));
        posy_min = bandWith+ry(ii); posy_max = ymax-(bandWith+ry(ii));
        if ii==1
            posx(ii) = posx_min+(posx_max-posx_min)*rand();
            posy(ii) = posy_min+(posy_max-posy_min)*rand();
        else
%             distEE=0;
%             p=0;
            maxDistEE = 0;
            for p=1:iMax
                posx(ii) = posx_min+(posx_max-posx_min)*rand();
                posy(ii) = posy_min+(posy_max-posy_min)*rand();
                [idx, distCC]= KDTreeBallQuery([posx(1:ii-1),posy(1:ii-1)],[posx(ii),posy(ii)],distMin+5*rmax*2);
%                 KDTreeClosestPoint([posx(1:ii-1),posy(1:ii-1)],[posx(ii),posy(ii)]); %%distCC=center-to-center spacing
                %angle
                idx = idx{1};
                distCC = distCC{1};
                if isempty(idx)
                    break
                else
                    angle = atan2(posy(ii)-posy(idx),posx(ii)-posx(idx));
                    % in elipse, (x,y) = (a cos(theta), b sin(theta)) where a
                    % and be are rs and ry
                    distEE = distCC-((rx(idx)'.*cos(angle)).^2+(ry(idx)'.*sin(angle)).^2).^0.5...
                        - ((rx(ii).*cos(angle)).^2+(ry(ii).*sin(angle)).^2).^0.5;
                    [minDistEE,~]=min(distEE);
                    if minDistEE>distMin
                        break
                    else
                        if minDistEE>maxDistEE
                            oldPosx = posx(ii);
                            oldPosy = posy(ii);
                            maxDistEE = minDistEE;
                        end
                    end
                end
            end
            if p==iMax
                disp(['The ' num2str(ii) 'th adhesion could not be positioned by ' num2str(distMin) ' pixel apart. ' num2str(maxDistEE) ' close. Skipping this adhesion ...'])
                posDiscard = [posDiscard; ii];
%                 posx(ii) = oldPosx;
%                 posy(ii) = oldPosy;
            end
        end
    end
    %% displacement field
    ux = zeros(size(x_mat_u));
    uy = zeros(size(x_mat_u));
    force_x = zeros(size(x_mat_u));
    force_y = zeros(size(x_mat_u));
    % filtering misplaced adhesions
   n = n - length(posDiscard);
   posx(posDiscard) = [];
   posy(posDiscard) = [];
   rx(posDiscard) = [];
   ry(posDiscard) = [];
   f(posDiscard) = [];
    
    for ii=1:n
        [ux_cur, uy_cur]= fwdSolution(x_mat_u,y_mat_u,E,xmin,xmax,ymin,ymax,...
        @(x,y) assumedForceAniso2D(1,x,y,posx(ii),posy(ii),0,f(ii),rx(ii),ry(ii),forceType),...
        @(x,y) assumedForceAniso2D(2,x,y,posx(ii),posy(ii),0,f(ii),rx(ii),ry(ii),forceType),'fft',[],meshPtsFwdSol);
        ux = ux + ux_cur;
        uy = uy + uy_cur;

        force_x_cur = assumedForceAniso2D(1,x_mat_u,y_mat_u,posx(ii),posy(ii),0,f(ii),rx(ii),ry(ii),forceType);
        force_y_cur = assumedForceAniso2D(2,x_mat_u,y_mat_u,posx(ii),posy(ii),0,f(ii),rx(ii),ry(ii),forceType);
        force_x = force_x + force_x_cur;
        force_y = force_y + force_y_cur;
    end
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
    MD.notes_=['Created for single force test purposes with f=' num2str(fmax) ' and n=' num2str(n)]; 

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
    params.maxFlowSpeed = 20;
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
        params.regParam = 2.4e-4;
    elseif strcmp(method,'L2')
        params.solMethodBEM = 'QR';
        params.regParam = 1e-7;
    else
        display('The method should be either L1 or L2. The input does not belong to any of those. We use L2 as a default.')
        params.solMethodBEM = 'QR';
    end
    params.method = 'FastBEM';
    params.useLcurve = true;
    params.basisClassTblPath = '/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM basis functions/basisClass8kPaSimul.mat';
    MD.getPackage(iPack).getProcess(4).setPara(params);
    MD.getPackage(iPack).getProcess(4).run();

    MD.save;
forceField=MD.getPackage(iPack).getProcess(4).loadChannelOutput;
else
    orgPath=[baseDataPath filesep 'Original'];
    analysisFolder = baseDataPath;
    if ~exist(dataPath,'dir')
        mkdir(dataPath);
    end
    %% loading original displacement field and force field
    orgMethod = method;
    backupDP = dataPath;
    load([orgPath filesep 'data.mat']);
    method = orgMethod;
    dataPath = backupDP;
    %% Now force reconstruction via movieData (non-GUI mode)
    % Retrieve current location
    MD=MovieData.load(fullfile(analysisFolder,'movieData.mat'));

    %% Create TFM package and retrieve package index
    iPack=  MD.getPackageIndex('TFMPackage');
    %% Create force reconstruction process and run
    params = MD.getPackage(iPack).getProcess(4).funParams_;
    if strcmp(method,'L1')
        params.solMethodBEM = '1NormReg';
        params.regParam = 3.8e-5;
    elseif strcmp(method,'L2')
        params.solMethodBEM = 'QR';
        params.regParam = 4.4e-7;
    else
        display('The method should be either L1 or L2. The input does not belong to any of those. We use L2 as a default.')
        params.solMethodBEM = 'QR';
    end
    MD.getPackage(iPack).getProcess(4).setPara(params);

    MD.getPackage(iPack).getProcess(4).run();
    forceField=MD.getPackage(iPack).getProcess(4).loadChannelOutput;
    mkdir([dataPath filesep 'TFMPackage/forceField']);
    save([dataPath filesep 'TFMPackage/forceField/forceField.mat'], 'forceField')
    save([dataPath filesep 'movieData.mat'], 'MD');
end

%% Postprocessing - saving and analyzing force field
% Loading displacement field and force field
% Load the displField
disp('Calculating RMS error and error on adhesions (EOAs) ... ')
% Load the forcefield
[fMap,XI,YI]=generateHeatmapFromField(forceField);
cropInfo = [XI(1,1), YI(1,1)];
% RMS error
% get x-component and y-component of forceField
beta = [forceField.vec(:,1); forceField.vec(:,2)];
load(MD.getPackage(iPack).getProcess(4).outFilePaths_{2},'forceMesh')
[fx,fy]=calcForcesFromCoef(forceMesh,beta,XI,YI,'new');
% compare fx, fy with force_x and force_y 
force_x_crop = force_x(cropInfo(2):cropInfo(2)+size(YI,1)-1,cropInfo(1):cropInfo(1)+size(XI,2)-1);
force_y_crop = force_y(cropInfo(2):cropInfo(2)+size(YI,1)-1,cropInfo(1):cropInfo(1)+size(XI,2)-1);

rmsError = ((force_x_crop-fx).^2-(force_y_crop-fy).^2).^0.5;
rmsError = sum(rmsError(:))/(size(rmsError,1)*size(rmsError,2));

% EOA: error on adhesion
% filter making with adhesion location
mask = false(size(XI));
distThres = 4;
for ii=1:n
    mask = mask | ((XI-posx(ii)).^2+(YI-posy(ii)).^2).^0.5<=distThres;
end
% EOA calculation
EOA = ((force_x_crop.*mask-fx.*mask).^2-(force_y_crop.*mask-fy.*mask).^2).^0.5;
EOA = sum(EOA(:))/floor(pi*distThres^2/4)/n; % error per adhesion
%% L-curve scrutinization - regularization parameter, convexness?
lcurvePath = MD.getPackage(iPack).getProcess(4).outFilePaths_{3};
%% save original forcefield
[h2,orgMap] = generateHeatmapFromGridData(x_mat_u,y_mat_u,force_x,force_y,[dataPath '/Original forcefield'],0,fmax,false,size(XI,2),size(YI,1));
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
